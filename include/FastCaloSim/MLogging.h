/*
  Copyright (C) 2002-2024 CERN for the benefit of the ATLAS collaboration
*/

#ifndef ISF_FASTCALOSIMEVENT_MLogging_h
#define ISF_FASTCALOSIMEVENT_MLogging_h

#include <FastCaloSim/FastCaloSim_export.h>

#include <iomanip>
#include <iostream>

#include <TNamed.h>  //for ClassDef
namespace MSG
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
}  // end namespace MSG

// Macro for use outside classes.
// Use this in standalone functions or static methods.
#define ATH_MSG_NOCLASS(logger_name, x) \
  do { \
    logger_name.msg() << logger_name.startMsg(MSG::ALWAYS, __FILE__, __LINE__) \
                      << x << std::endl; \
  } while (0)

// Declare the class accessories in a namespace
// namespace ISF_FCS {
//__attribute__((unused)) static const char *LevelNames[MSG::NUM_LEVELS] = {
//    "NIL", "VERBOSE", "DEBUG", "INFO", "WARNING", "ERROR", "FATAL", "ALWAYS"};
//} // end namespace ISF_FCS

// Declare the class in a namespace
namespace ISF_FCS
{

// We can define a number of macros here to replace standard ATH_MSG
// macros. This can only be done outside Athena or the compiler complains.
typedef std::ostream MsgStream;

#define ATH_MSG_LVL(enum_lvl, x) \
  do { \
    if (this->msgLvl(enum_lvl)) \
      this->msg() << this->startMsg(enum_lvl, __FILE__, __LINE__) << x \
                  << std::endl; \
  } while (0)

#define ATH_MSG_LVL_NOCHK(enum_lvl, x) \
  do { \
    this->msg() << this->startMsg(enum_lvl, __FILE__, __LINE__) << x \
                << std::endl; \
  } while (0)

#define ATH_MSG_VERBOSE(x) ATH_MSG_LVL(MSG::VERBOSE, x)
#define ATH_MSG_DEBUG(x) ATH_MSG_LVL(MSG::DEBUG, x)
#define ATH_MSG_INFO(x) ATH_MSG_LVL_NOCHK(MSG::INFO, x)
#define ATH_MSG_WARNING(x) ATH_MSG_LVL_NOCHK(MSG::WARNING, x)
#define ATH_MSG_ERROR(x) ATH_MSG_LVL_NOCHK(MSG::ERROR, x)
#define ATH_MSG_FATAL(x) ATH_MSG_LVL_NOCHK(MSG::FATAL, x)

// Set up a stream that can be used like: ATH_MSG(INFO) << "hello" <<
// END_MSG(INFO); It needs to only write the left columns once, until it is fed
// another END_MSG

// Provide a stream
#define ATH_MSG(lvl) this->stream(MSG::lvl, __FILE__, __LINE__)
// Add a new line if the level is in use, and end any stream
#define END_MSG(lvl) this->streamerEndLine(MSG::lvl)
// Force a new line, and end any stream
#define endmsg this->streamerEndLine(MSG::INFO)

class FASTCALOSIM_EXPORT MLogging
{
public:
  /// Constructor
  MLogging() {};
  /// Copy constructor
  MLogging(const MLogging& other)
      : m_level(other.m_level) {};
  /// Assignment operator
  MLogging& operator=(MLogging other)
  {
    m_level = other.m_level;
    return *this;
  };

  /// Destructor
  virtual ~MLogging() {};

  /// Retrieve output level
  MSG::Level level() const { return m_level; }
  /// Update outputlevel
  virtual void setLevel(int level);

  /// Make a message to decorate the start of logging
  static std::string startMsg(MSG::Level lvl,
                              const std::string& file,
                              int line);

  /// Return a stream for sending messages directly (no decoration)
  MsgStream& msg() const { return *m_msg; }
  /// Return a stream for sending messages (incomplete decoration)
  MsgStream& msg(const MSG::Level lvl) const;
  /// Return a decorated starting stream for sending messages
  MsgStream& stream(MSG::Level lvl, std::string file, int line) const;
  /// Check whether the logging system is active at the provided verbosity level
  bool msgLvl(const MSG::Level lvl) const;

  /// Print a whole decorated log message and then end the line.
  void print(MSG::Level lvl,
             std::string file,
             int line,
             std::string message) const;

  /// Update and end the line if we print this level
  std::string streamerEndLine(MSG::Level lvl) const;

private:
  /// Checking the state of the streamer.
  bool streamerInLine() const { return m_streamer_in_line; }
  /// Update if a new start is happening.
  void streamerInLine(bool is_in_line) const;
  /// Check if a new start should be done (changed file or level)
  bool streamerNeedStart(MSG::Level lvl, std::string file) const;

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

}  // namespace ISF_FCS

#endif  // End header guard
