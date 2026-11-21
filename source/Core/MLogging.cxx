// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#include <array>
#include <atomic>
#include <format>
#include <iostream>
#include <memory>
#include <mutex>
#include <sstream>
#include <string>
#include <string_view>
#include <utility>

#include "FastCaloSim/Core/MLogging.h"

namespace
{
// MLogging was called ISF_FCS::MLogging before the library was split off from
// Athena, and parametrization files written back then reference that name in
// their streamer info. Registering the old name as an alternate lets ROOT
// resolve it to the renamed class so legacy files keep loading without
// warnings.
const ROOT::TClassAlt* dummy =
    ROOT::AddClassAlternate("FastCaloSim::MLogging", "ISF_FCS::MLogging");

// Writes each message as one decorated line to std::cout. The mutex keeps
// lines from different threads apart.
class CoutSink final : public FastCaloSim::LogSink
{
public:
  void write(FastCaloSim::MSG::Level lvl,
             std::string_view file,
             int line,
             std::string_view text) override
  {
    const std::string out = std::format(
        "{}{}\n", FastCaloSim::MLogging::startMsg(lvl, file, line), text);
    const std::scoped_lock lock(m_mutex);
    std::cout << out << std::flush;
  }

private:
  std::mutex m_mutex;
};

// The default sink and the slot of the current sink are never destroyed, so
// that messages from destructors at exit still work.
auto default_sink() -> const std::shared_ptr<FastCaloSim::LogSink>&
{
  static const auto* sink =
      new std::shared_ptr<FastCaloSim::LogSink>(std::make_shared<CoutSink>());
  return *sink;
}

auto sink_slot() -> std::atomic<std::shared_ptr<FastCaloSim::LogSink>>&
{
  static auto* slot =
      new std::atomic<std::shared_ptr<FastCaloSim::LogSink>>(default_sink());
  return *slot;
}

// Called from destructors, so it must not throw. A message that the installed
// sink fails to take goes to the default sink instead, so that no message is
// lost. Newlines at the end of the text are removed: sinks end lines
// themselves.
void send(FastCaloSim::MSG::Level lvl,
          std::string_view file,
          int line,
          std::string_view text) noexcept
{
  while (text.ends_with('\n')) {
    text.remove_suffix(1);
  }
  try {
    sink_slot().load()->write(lvl, file, line, text);
    return;
  } catch (...) {
  }
  try {
    default_sink()->write(lvl, file, line, text);
  } catch (...) {
  }
}

// The line that the FCS_MSG(lvl) << ... statements of one thread build until
// END_FCS_MSG.
struct ThreadState
{
  bool open = false;
  const FastCaloSim::MLogging* owner = nullptr;
  FastCaloSim::MSG::Level lvl = FastCaloSim::MSG::NIL;
  std::string_view file;
  int line = 0;
  std::ostringstream text;

  void end_line()
  {
    if (!open) {
      return;
    }
    open = false;
    send(lvl, file, line, text.view());
    text.str({});
  }
};

// Plain thread_local values have no destructor and stay usable while the
// thread exits, unlike a thread_local ThreadState: static destructors that log
// run after the thread_local objects of the main thread are gone.
thread_local ThreadState* t_state = nullptr;
thread_local bool t_exiting = false;

struct ThreadStateOwner
{
  ThreadStateOwner() = default;
  ThreadStateOwner(const ThreadStateOwner&) = delete;
  auto operator=(const ThreadStateOwner&) -> ThreadStateOwner& = delete;
  ~ThreadStateOwner()
  {
    t_exiting = true;
    if (ThreadState* state = std::exchange(t_state, nullptr)) {
      state->end_line();
      delete state;
    }
  }
};

auto thread_state() -> ThreadState&
{
  if (t_state == nullptr) {
    t_state = new ThreadState;
    // Deleted when the thread exits. State created while the thread exits
    // (logging from later destructors) is not deleted.
    if (!t_exiting) {
      [[maybe_unused]] thread_local ThreadStateOwner owner;
    }
  }
  return *t_state;
}
}  // namespace

namespace FastCaloSim
{

void set_log_sink(std::shared_ptr<LogSink> sink)
{
  sink_slot().store(sink ? std::move(sink) : default_sink());
}

LogMessage::LogMessage(MSG::Level lvl, const char* file, int line)
    : m_lvl(lvl)
    , m_file(file)
    , m_line(line)
{
}

LogMessage::~LogMessage()
{
  send(m_lvl, m_file, m_line, m_text.view());
}

/// Update outputlevel
void MLogging::setLevel(int level)
{
  level = (level >= MSG::NUM_LEVELS) ? MSG::ALWAYS
      : (level < MSG::NIL)           ? MSG::NIL
                                     : level;
  m_level = MSG::Level(level);
}

LineStream::LineStream(const MLogging* owner,
                       MSG::Level lvl,
                       const char* file,
                       int line,
                       bool active)
    : m_owner(owner)
    , m_lvl(lvl)
    , m_file(file)
    , m_line(line)
{
  if (active) {
    m_text.emplace();
  }
}

/// Add the text of the statement to the open line of this thread. Another
/// object, level or file ends the line and starts a new one.
LineStream::~LineStream()
{
  if (!m_text) {
    return;
  }
  ThreadState& state = thread_state();
  if (state.open
      && (m_owner != state.owner || m_lvl != state.lvl
          || std::string_view(m_file) != state.file))
  {
    state.end_line();
  }
  if (!state.open) {
    state.open = true;
    state.owner = m_owner;
    state.lvl = m_lvl;
    state.file = m_file;
    state.line = m_line;
  }
  state.text << m_text->view();
  // While the thread exits, nothing would end the line later.
  if (m_end || t_exiting) {
    state.end_line();
  }
}

auto MLogging::stream(MSG::Level lvl, const char* file, int line) const
    -> LineStream
{
  return {this, lvl, file, line, this->msgLvl(lvl)};
}

auto MLogging::streamerEndLine(MSG::Level lvl) const -> LineStream::End
{
  return {this->msgLvl(lvl)};
}

/// Print a message for the start of logging
auto MLogging::startMsg(MSG::Level lvl, std::string_view file, int line)
    -> std::string
{
  static constexpr std::array<const char*, MSG::NUM_LEVELS> level_names = {
      "NIL", "VERBOSE", "DEBUG", "INFO", "WARNING", "ERROR", "FATAL", "ALWAYS"};
  // File name from the last '/', or its last 20 characters if longer
  constexpr std::size_t name_len = 20;
  const auto last_slash = file.find_last_of('/');
  std::size_t trim_point =
      last_slash == std::string_view::npos ? 0 : last_slash;
  if (file.length() - trim_point > name_len) {
    trim_point = file.length() - name_len;
  }
  return std::format("{:>20}:{:<5}{:>10}",
                     file.substr(trim_point),
                     line,
                     std::format("({}) ", level_names[lvl]));
}

}  // namespace FastCaloSim
