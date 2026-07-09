// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#pragma once

#include <iomanip>
#include <iostream>

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

// Macro for use outside classes.
// Use this in standalone functions or static methods.
#define FCS_MSG_NOCLASS(logger_name, x) \
  do { \
    logger_name.msg() << logger_name.startMsg( \
        ::FastCaloSim::MSG::ALWAYS, __FILE__, __LINE__) \
                      << x << std::endl; \
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
      this->msg() << this->startMsg(enum_lvl, __FILE__, __LINE__) << x \
                  << std::endl; \
  } while (0)

#define FCS_MSG_LVL_NOCHK(enum_lvl, x) \
  do { \
    this->msg() << this->startMsg(enum_lvl, __FILE__, __LINE__) << x \
                << std::endl; \
  } while (0)

#define FCS_MSG_VERBOSE(x) FCS_MSG_LVL(::FastCaloSim::MSG::VERBOSE, x)
#define FCS_MSG_DEBUG(x) FCS_MSG_LVL(::FastCaloSim::MSG::DEBUG, x)
#define FCS_MSG_INFO(x) FCS_MSG_LVL_NOCHK(::FastCaloSim::MSG::INFO, x)
#define FCS_MSG_WARNING(x) FCS_MSG_LVL_NOCHK(::FastCaloSim::MSG::WARNING, x)
#define FCS_MSG_ERROR(x) FCS_MSG_LVL_NOCHK(::FastCaloSim::MSG::ERROR, x)
#define FCS_MSG_FATAL(x) FCS_MSG_LVL_NOCHK(::FastCaloSim::MSG::FATAL, x)

// Set up a stream that can be used like: FCS_MSG(INFO) << "hello" <<
// END_FCS_MSG(INFO); It needs to only write the left columns once, until it is
// fed another END_FCS_MSG

// Provide a stream
#define FCS_MSG(lvl) this->stream(::FastCaloSim::MSG::lvl, __FILE__, __LINE__)
// Add a new line if the level is in use, and end any stream
#define END_FCS_MSG(lvl) this->streamerEndLine(::FastCaloSim::MSG::lvl)
// Force a new line, and end any stream
#define FCS_ENDMSG this->streamerEndLine(::FastCaloSim::MSG::INFO)

class FASTCALOSIM_EXPORT MLogging
{
public:
  /// Constructor
  MLogging() = default;
  /// Copy constructor
  MLogging(const MLogging& other)
      : m_level(other.m_level) {};
  /// Assignment operator
  auto operator=(MLogging other) -> MLogging&
  {
    m_level = other.m_level;
    return *this;
  };

  /// Destructor
  virtual ~MLogging() = default;

  /// Retrieve output level
  auto level() const -> MSG::Level { return m_level; }
  /// Update outputlevel
  virtual void setLevel(int level);

  /// Make a message to decorate the start of logging
  static auto startMsg(MSG::Level lvl, const std::string& file, int line)
      -> std::string;

  /// Return a stream for sending messages directly (no decoration)
  auto msg() const -> MsgStream& { return *m_msg; }
  /// Return a stream for sending messages (incomplete decoration)
  auto msg(const MSG::Level lvl) const -> MsgStream&;
  /// Return a decorated starting stream for sending messages
  auto stream(MSG::Level lvl, std::string file, int line) const -> MsgStream&;
  /// Check whether the logging system is active at the provided verbosity level
  auto msgLvl(const MSG::Level lvl) const -> bool;

  /// Print a whole decorated log message and then end the line.
  void print(MSG::Level lvl,
             std::string file,
             int line,
             std::string message) const;

  /// Update and end the line if we print this level
  auto streamerEndLine(MSG::Level lvl) const -> std::string;

private:
  /// Checking the state of the streamer.
  auto streamerInLine() const -> bool { return m_streamer_in_line; }
  /// Update if a new start is happening.
  void streamerInLine(bool is_in_line) const;
  /// Check if a new start should be done (changed file or level)
  auto streamerNeedStart(MSG::Level lvl, std::string file) const -> bool;

  MSG::Level m_level = MSG::INFO;  //! Do not persistify!

  MsgStream* m_msg = &std::cout;  //! Do not persistify!
  MsgStream m_null_msg = MsgStream(nullptr);  //! Do not persistify!
  MsgStream* m_null_msg_ptr = &m_null_msg;  //! Do not persistify!

  mutable bool m_streamer_in_line = false;  //! Do not persistify!
  mutable MSG::Level m_streamer_has_lvl = MSG::NIL;  //! Do not persistify!
  mutable std::string m_streamer_from_file = "";  //! Do not persistify!

  // Version number 0 to tell ROOT not to store this.
  ClassDef(MLogging, 0)
};
}  // namespace FastCaloSim
