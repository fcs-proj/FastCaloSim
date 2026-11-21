// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#pragma once

#include <memory>
#include <optional>
#include <ostream>
#include <sstream>
#include <string>
#include <string_view>

#include <FastCaloSim/FastCaloSim_export.h>
#include <TNamed.h>  //for ClassDef

namespace FastCaloSim::MSG
{
enum Level
{
  NIL = 0,
  VERBOSE,
  DEBUG,
  INFO,
  WARNING,
  ERROR,
  FATAL,
  ALWAYS,
  NUM_LEVELS
};  // enum Level
}  // namespace FastCaloSim::MSG

namespace FastCaloSim
{
/// Receives every log message that passes the level check. Install one with
/// set_log_sink() to forward FastCaloSim messages to another logging system,
/// for example a framework's message service.
class FASTCALOSIM_EXPORT LogSink
{
public:
  virtual ~LogSink() = default;
  /// One message: `file` and `line` locate the call site, `text` is the
  /// message as written at the call site, without decoration and without an
  /// added newline. Called concurrently from all threads that log. An exception
  /// thrown here does not propagate: the message then goes to the default sink.
  virtual void write(MSG::Level lvl,
                     std::string_view file,
                     int line,
                     std::string_view text) = 0;
};

/// Replace the log sink for all threads; nullptr restores the default sink,
/// which writes each message as one decorated line to std::cout.
FASTCALOSIM_EXPORT void set_log_sink(std::shared_ptr<LogSink> sink);

/// Collects the text of one message and passes it to the log sink when it goes
/// out of scope. The FCS_MSG_* macros create one per message.
class FASTCALOSIM_EXPORT LogMessage
{
public:
  LogMessage(MSG::Level lvl, const char* file, int line);
  LogMessage(const LogMessage&) = delete;
  auto operator=(const LogMessage&) -> LogMessage& = delete;
  ~LogMessage();

  auto stream() -> std::ostream& { return m_text; }

private:
  MSG::Level m_lvl;
  const char* m_file;
  int m_line;
  std::ostringstream m_text;
};

class MLogging;

/// Collects the text of one FCS_MSG(lvl) << ... statement and adds it to the
/// open line of this thread when the statement ends. A statement nested in
/// another one (a function called in the expression that streams itself)
/// therefore finishes first and does not mix with the outer text.
class FASTCALOSIM_EXPORT LineStream
{
public:
  /// Streamed by END_FCS_MSG: the line ends with this statement
  struct [[nodiscard]] End
  {
    bool active;
  };

  LineStream(const MLogging* owner,
             MSG::Level lvl,
             const char* file,
             int line,
             bool active);
  LineStream(const LineStream&) = delete;
  auto operator=(const LineStream&) -> LineStream& = delete;
  ~LineStream();

  template<typename T>
  auto operator<<(const T& value) -> LineStream&
  {
    if (m_text)
      *m_text << value;
    return *this;
  }
  auto operator<<(std::ostream& (*manipulator)(std::ostream&)) -> LineStream&
  {
    if (m_text)
      *m_text << manipulator;
    return *this;
  }
  auto operator<<(End end) -> LineStream&
  {
    m_end = m_end || end.active;
    return *this;
  }

private:
  const MLogging* m_owner;
  MSG::Level m_lvl;
  const char* m_file;
  int m_line;
  bool m_end = false;
  std::optional<std::ostringstream> m_text;  // only if the level is active
};
}  // namespace FastCaloSim

// Macro for use outside classes.
// Use this in standalone functions or static methods.
#define FCS_MSG_NOCLASS(logger_name, x) \
  do { \
    static_cast<void>(logger_name); \
    ::FastCaloSim::LogMessage fcs_log_message( \
        ::FastCaloSim::MSG::ALWAYS, __FILE__, __LINE__); \
    fcs_log_message.stream() << x; \
  } while (0)

namespace FastCaloSim
{
// We can define a number of macros here to replace the standard Athena
// messaging macros. This can only be done outside Athena or the compiler
// complains.
using MsgStream = std::ostream;

#define FCS_MSG_LVL(enum_lvl, x) \
  do { \
    if (this->msgLvl(enum_lvl)) \
      FCS_MSG_LVL_NOCHK(enum_lvl, x); \
  } while (0)

#define FCS_MSG_LVL_NOCHK(enum_lvl, x) \
  do { \
    ::FastCaloSim::LogMessage fcs_log_message(enum_lvl, __FILE__, __LINE__); \
    fcs_log_message.stream() << x; \
  } while (0)

#define FCS_MSG_VERBOSE(x) FCS_MSG_LVL(::FastCaloSim::MSG::VERBOSE, x)
#define FCS_MSG_DEBUG(x) FCS_MSG_LVL(::FastCaloSim::MSG::DEBUG, x)
#define FCS_MSG_INFO(x) FCS_MSG_LVL_NOCHK(::FastCaloSim::MSG::INFO, x)
#define FCS_MSG_WARNING(x) FCS_MSG_LVL_NOCHK(::FastCaloSim::MSG::WARNING, x)
#define FCS_MSG_ERROR(x) FCS_MSG_LVL_NOCHK(::FastCaloSim::MSG::ERROR, x)
#define FCS_MSG_FATAL(x) FCS_MSG_LVL_NOCHK(::FastCaloSim::MSG::FATAL, x)

// Build one line in pieces: FCS_MSG(INFO) << "a"; ... FCS_MSG(INFO) << "b"
// << END_FCS_MSG(INFO); Each thread has one open line. It is passed to the log
// sink at the end of the statement that streams END_FCS_MSG, or when the same
// thread streams from another object, at another level or from another file.
// FCS_MSG_* messages in between are passed on by themselves and leave the line
// open.

// Provide a stream for one statement
#define FCS_MSG(lvl) this->stream(::FastCaloSim::MSG::lvl, __FILE__, __LINE__)
// Stream this to end the line if the level is in use
#define END_FCS_MSG(lvl) this->streamerEndLine(::FastCaloSim::MSG::lvl)
// Stream this to end the line
#define FCS_ENDMSG this->streamerEndLine(::FastCaloSim::MSG::INFO)

class FASTCALOSIM_EXPORT MLogging
{
public:
  /// Constructor
  MLogging() = default;
  MLogging(const MLogging&) = default;
  auto operator=(const MLogging&) -> MLogging& = default;

  /// Destructor
  virtual ~MLogging() = default;

  /// Retrieve output level
  auto level() const -> MSG::Level { return m_level; }
  /// Update outputlevel
  virtual void setLevel(int level);

  /// Make a message to decorate the start of logging
  static auto startMsg(MSG::Level lvl, std::string_view file, int line)
      -> std::string;

  /// Return a stream for one statement that adds to the open line (FCS_MSG)
  auto stream(MSG::Level lvl, const char* file, int line) const -> LineStream;
  /// Check whether the logging system is active at the provided verbosity
  /// level: VERBOSE and DEBUG depend on the output level, all other levels are
  /// always printed
  auto msgLvl(const MSG::Level lvl) const -> bool
  {
    if (lvl == MSG::VERBOSE || lvl == MSG::DEBUG)
      return m_level <= lvl;
    return true;
  }

  /// Marker that ends the line if we print this level (END_FCS_MSG)
  auto streamerEndLine(MSG::Level lvl) const -> LineStream::End;

private:
  MSG::Level m_level = MSG::INFO;  //! Do not persistify!

  // Version number 0 to tell ROOT not to store this.
  ClassDef(MLogging, 0)
};
}  // namespace FastCaloSim
