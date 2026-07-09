// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

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
    return *m_null_msg_ptr;
  // This is different to the last stream.
  if (this->streamerNeedStart(lvl, file)) {
    // end the open stream
    if (this->streamerInLine())
      std::cout << this->streamerEndLine(lvl);
    this->msg() << this->startMsg(lvl, file, line);
    this->streamerInLine(true);
  }
  return *m_msg;
}

/// Check whether the logging system is active at the provided verbosity level
auto MLogging::msgLvl(const MSG::Level lvl) const -> bool
{
  if (lvl == MSG::VERBOSE || lvl == MSG::DEBUG)
    return m_level <= lvl;
  // All other messages print always
  return true;
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
    *m_msg << std::endl;
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
  int col1_len = 20;
  int col2_len = 5;
  int col3_len = 10;
  auto last_slash = file.find_last_of('/');
  int path_len = last_slash == std::string::npos ? 0 : last_slash;
  int trim_point = path_len;
  int total_len = file.length();
  if (total_len - path_len > col1_len)
    trim_point = total_len - col1_len;
  std::string trimmed_name = file.substr(trim_point);
  const char* LevelNames[MSG::NUM_LEVELS] = {
      "NIL", "VERBOSE", "DEBUG", "INFO", "WARNING", "ERROR", "FATAL", "ALWAYS"};
  std::string level = LevelNames[lvl];
  std::string level_string = std::string("(") + level + ") ";
  std::stringstream output;
  output << std::setw(col1_len) << std::right << trimmed_name << ":"
         << std::setw(col2_len) << std::left << line << std::setw(col3_len)
         << std::right << level_string;
  return output.str();
}

}  // namespace FastCaloSim
