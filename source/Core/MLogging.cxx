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

// Stream for messages below the output level: it has no buffer, so everything
// written to it is discarded. One per thread, as writing sets its state flags.
auto null_stream() -> std::ostream&
{
  thread_local std::ostream stream(nullptr);
  return stream;
}

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

// Never destroyed, so that messages from destructors at exit still work.
auto sink_slot() -> std::atomic<std::shared_ptr<FastCaloSim::LogSink>>&
{
  static auto* slot = new std::atomic<std::shared_ptr<FastCaloSim::LogSink>>(
      std::make_shared<CoutSink>());
  return *slot;
}

void send(FastCaloSim::MSG::Level lvl,
          std::string_view file,
          int line,
          std::string_view text)
{
  sink_slot().load()->write(lvl, file, line, text);
}

// The line of this thread that FCS_MSG(lvl) << ... builds until END_FCS_MSG.
struct OpenLine
{
  bool open = false;
  FastCaloSim::MSG::Level lvl = FastCaloSim::MSG::NIL;
  std::string_view file;
  int line = 0;
  std::ostringstream text;

  void end()
  {
    if (!open)
      return;
    open = false;
    send(lvl, file, line, text.view());
  }
  ~OpenLine() { end(); }
};

auto open_line() -> OpenLine&
{
  thread_local OpenLine line;
  return line;
}
}  // namespace

namespace FastCaloSim
{

void set_log_sink(std::shared_ptr<LogSink> sink)
{
  if (!sink)
    sink = std::make_shared<CoutSink>();
  sink_slot().store(std::move(sink));
}

LogMessage::LogMessage(MSG::Level lvl, const char* file, int line)
    : m_lvl(lvl)
    , m_file(file)
    , m_line(line)
{
}

LogMessage::~LogMessage()
{
  // A message ends any line this thread was building, which comes first.
  open_line().end();
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

/// Return the open line of this thread. A different level or file ends the
/// line and starts a new one.
auto MLogging::stream(MSG::Level lvl, const char* file, int line) const
    -> MsgStream&
{
  // If we shouldn't print this just return a dummy stream.
  if (!this->msgLvl(lvl))
    return null_stream();
  OpenLine& open = open_line();
  if (open.open && (lvl != open.lvl || file != open.file))
    open.end();
  if (!open.open) {
    open.open = true;
    open.lvl = lvl;
    open.file = file;
    open.line = line;
    open.text.str({});
  }
  return open.text;
}

/// End the open line if we print this level
auto MLogging::streamerEndLine(MSG::Level lvl) const -> std::string
{
  if (this->msgLvl(lvl))
    open_line().end();
  return "";
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
  if (file.length() - trim_point > name_len)
    trim_point = file.length() - name_len;
  return std::format("{:>20}:{:<5}{:>10}",
                     file.substr(trim_point),
                     line,
                     std::format("({}) ", level_names[lvl]));
}

}  // namespace FastCaloSim
