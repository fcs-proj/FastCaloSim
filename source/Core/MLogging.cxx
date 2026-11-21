// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#include <array>
#include <format>

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
}  // namespace

namespace FastCaloSim
{

/// Update outputlevel
void MLogging::setLevel(int level)
{
  level = (level >= MSG::NUM_LEVELS) ? MSG::ALWAYS
      : (level < MSG::NIL)           ? MSG::NIL
                                     : level;
  m_level = MSG::Level(level);
}

// startMsg defined at base of file.

/// Return a stream for sending messages (incomplete decoration)
auto MLogging::msg(const MSG::Level lvl) const -> MsgStream&
{
  return this->stream(lvl, "", -1);
};

/// Open a stream for sending log messages
// check the lvl of logging is active, otherwise provide dummy stream.
// if the proposed streamer doesn't match the current one
// end any running lines and start a new decorated line
// provide the stream at the end
auto MLogging::stream(MSG::Level lvl, std::string file, int line) const
    -> MsgStream&
{
  // If we shouldn't print this just return a dummy stream.
  if (!this->msgLvl(lvl))
    return null_stream();
  // This is different to the last stream.
  if (this->streamerNeedStart(lvl, file)) {
    // end the open stream
    if (this->streamerInLine())
      std::cout << this->streamerEndLine(lvl);
    this->msg() << this->startMsg(lvl, file, line);
    this->streamerInLine(true);
  }
  return this->msg();
}

/// Print a whole log message and then end the line.
void MLogging::print(MSG::Level lvl,
                     std::string file,
                     int line,
                     std::string message) const
{
  if (this->msgLvl(lvl)) {
    this->stream(lvl, file, line) << message << this->streamerEndLine(lvl);
  }
}

/// Update and end the line if we print this level
auto MLogging::streamerEndLine(MSG::Level lvl) const -> std::string
{
  if (this->msgLvl(lvl)) {
    m_streamer_in_line = false;
    this->msg() << std::endl;
  }
  return "";
}

/// Update if a new start is happening.
void MLogging::streamerInLine(bool is_in_line) const
{
  m_streamer_in_line = is_in_line;
}

/// Check if a new start should be done (changed file or level)
auto MLogging::streamerNeedStart(MSG::Level lvl, std::string file) const -> bool
{
  // Are we in the middle of a stream of the same level from the same file.
  if (lvl == m_streamer_has_lvl && file == m_streamer_from_file
      && m_streamer_in_line)
    return false;
  // Otherwise time for a new start.
  m_streamer_from_file = file;
  m_streamer_has_lvl = lvl;
  return true;
}

// This is the same either way.
/// Print a message for the start of logging
auto MLogging::startMsg(MSG::Level lvl, const std::string& file, int line)
    -> std::string
{
  static constexpr std::array<const char*, MSG::NUM_LEVELS> level_names = {
      "NIL", "VERBOSE", "DEBUG", "INFO", "WARNING", "ERROR", "FATAL", "ALWAYS"};
  // File name from the last '/', or its last 20 characters if longer
  constexpr std::size_t name_len = 20;
  const auto last_slash = file.find_last_of('/');
  std::size_t trim_point = last_slash == std::string::npos ? 0 : last_slash;
  if (file.length() - trim_point > name_len)
    trim_point = file.length() - name_len;
  return std::format("{:>20}:{:<5}{:>10}",
                     file.substr(trim_point),
                     line,
                     std::format("({}) ", level_names[lvl]));
}

}  // namespace FastCaloSim
