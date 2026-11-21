// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <string>
#include <string_view>
#include <thread>
#include <vector>

#include <gtest/gtest.h>

#include "FastCaloSim/Core/MLogging.h"

namespace
{
using FastCaloSim::MSG::Level;

struct Record
{
  Level lvl;
  std::string file;
  int line;
  std::string text;
};

// Keeps every message it receives
class CaptureSink : public FastCaloSim::LogSink
{
public:
  void write(Level lvl,
             std::string_view file,
             int line,
             std::string_view text) override
  {
    const std::scoped_lock lock(m_mutex);
    m_records.push_back({lvl, std::string(file), line, std::string(text)});
  }

  auto records() -> std::vector<Record>
  {
    const std::scoped_lock lock(m_mutex);
    return m_records;
  }

private:
  std::mutex m_mutex;
  std::vector<Record> m_records;
};

// Fails on every message
class ThrowingSink : public FastCaloSim::LogSink
{
public:
  void write(Level /*lvl*/,
             std::string_view /*file*/,
             int /*line*/,
             std::string_view /*text*/) override
  {
    throw std::runtime_error("sink failure");
  }
};

// Writes the bare text to std::cerr, where death tests can match it
class CerrSink : public FastCaloSim::LogSink
{
public:
  void write(Level /*lvl*/,
             std::string_view /*file*/,
             int /*line*/,
             std::string_view text) override
  {
    std::cerr << text << '\n';
  }
};

// Logs through the macros like the FastCaloSim classes do
class Logger : public FastCaloSim::MLogging
{
public:
  void info(int i) const { FCS_MSG_INFO("info " << i); }
  void text(const char* text) const { FCS_MSG_INFO(text); }
  void debug(int i) const { FCS_MSG_DEBUG("debug " << i); }
  void warning(int i) const { FCS_MSG_WARNING("warning " << i); }
  void error(int i) const { FCS_MSG_ERROR("error " << i); }
  template<typename F>
  void debug_with(F f) const
  {
    FCS_MSG_DEBUG("debug " << f());
  }
  void streamed(int n) const
  {
    FCS_MSG(INFO) << "items:";
    for (int i = 0; i < n; ++i)
      FCS_MSG(INFO) << " " << i;
    FCS_MSG(INFO) << END_FCS_MSG(INFO);
  }
  void part(const char* text) const { FCS_MSG(INFO) << text; }
  void end() const { FCS_MSG(INFO) << END_FCS_MSG(INFO); }
  auto inner() const -> int
  {
    FCS_MSG_INFO("inner");
    return 1;
  }
  void nested() const
  {
    FCS_MSG(INFO) << "before " << inner() << " after" << END_FCS_MSG(INFO);
  }
  auto inner_line() const -> int
  {
    FCS_MSG(INFO) << "inner" << END_FCS_MSG(INFO);
    return 42;
  }
  void nested_line(const Logger& other) const
  {
    FCS_MSG(INFO) << "before " << other.inner_line() << " after"
                  << END_FCS_MSG(INFO);
  }
};

// Streams a line from its destructor, ended or not
class ExitLogger : public Logger
{
public:
  explicit ExitLogger(bool end_line)
      : m_end_line(end_line)
  {
  }
  ExitLogger(const ExitLogger&) = delete;
  auto operator=(const ExitLogger&) -> ExitLogger& = delete;
  ~ExitLogger() override
  {
    part("streamed at exit");
    if (m_end_line)
      end();
  }

private:
  bool m_end_line;
};

class LoggingTest : public ::testing::Test
{
protected:
  void SetUp() override { FastCaloSim::set_log_sink(m_sink); }
  void TearDown() override { FastCaloSim::set_log_sink(nullptr); }

  std::shared_ptr<CaptureSink> m_sink = std::make_shared<CaptureSink>();
};
}  // namespace

TEST_F(LoggingTest, MessageReachesSinkOnce)
{
  const Logger log;
  log.info(7);
  const auto records = m_sink->records();
  ASSERT_EQ(records.size(), 1U);
  EXPECT_EQ(records[0].lvl, FastCaloSim::MSG::INFO);
  EXPECT_EQ(records[0].text, "info 7");
  EXPECT_TRUE(records[0].file.ends_with("LoggingTests.cxx"));
  EXPECT_GT(records[0].line, 0);
}

// Every level follows the output level: at WARNING, INFO is not printed
TEST_F(LoggingTest, LevelFiltersAllLevels)
{
  Logger log;
  log.setLevel(FastCaloSim::MSG::WARNING);
  log.info(1);
  log.streamed(1);
  log.warning(2);
  log.error(3);
  const auto records = m_sink->records();
  ASSERT_EQ(records.size(), 2U);
  EXPECT_EQ(records[0].text, "warning 2");
  EXPECT_EQ(records[1].text, "error 3");
}

// Newlines at the end are removed, newlines inside are kept
TEST_F(LoggingTest, TrailingNewlinesRemoved)
{
  const Logger log;
  log.text("a\nb\n\n");
  log.text("\n");
  const auto records = m_sink->records();
  ASSERT_EQ(records.size(), 2U);
  EXPECT_EQ(records[0].text, "a\nb");
  EXPECT_EQ(records[1].text, "");
}

TEST_F(LoggingTest, LevelFiltersDebug)
{
  Logger log;
  log.debug(1);
  EXPECT_TRUE(m_sink->records().empty());
  log.setLevel(FastCaloSim::MSG::DEBUG);
  log.debug(2);
  const auto records = m_sink->records();
  ASSERT_EQ(records.size(), 1U);
  EXPECT_EQ(records[0].lvl, FastCaloSim::MSG::DEBUG);
  EXPECT_EQ(records[0].text, "debug 2");
}

TEST_F(LoggingTest, StreamedLineIsOneMessage)
{
  const Logger log;
  log.streamed(3);
  const auto records = m_sink->records();
  ASSERT_EQ(records.size(), 1U);
  EXPECT_EQ(records[0].text, "items: 0 1 2");
}

TEST_F(LoggingTest, DisabledMessageIsNotEvaluated)
{
  const Logger log;
  int calls = 0;
  log.debug_with([&calls] { return ++calls; });
  EXPECT_EQ(calls, 0);
  EXPECT_TRUE(m_sink->records().empty());
}

// A message in between is passed on by itself and leaves the line open
TEST_F(LoggingTest, MessageLeavesOpenLine)
{
  const Logger log;
  log.part("open");
  log.info(1);
  log.end();
  const auto records = m_sink->records();
  ASSERT_EQ(records.size(), 2U);
  EXPECT_EQ(records[0].text, "info 1");
  EXPECT_EQ(records[1].text, "open");
}

// Also when the message is logged while the line is being written
TEST_F(LoggingTest, NestedMessageKeepsStreamedText)
{
  const Logger log;
  log.nested();
  const auto records = m_sink->records();
  ASSERT_EQ(records.size(), 2U);
  EXPECT_EQ(records[0].text, "inner");
  EXPECT_EQ(records[1].text, "before 1 after");
}

// A line streamed in the expression of another one finishes first, from the
// same object or from another one
TEST_F(LoggingTest, NestedLinesKeepStreamedText)
{
  const Logger a;
  const Logger b;
  a.nested_line(a);
  a.nested_line(b);
  const auto records = m_sink->records();
  ASSERT_EQ(records.size(), 4U);
  EXPECT_EQ(records[0].text, "inner");
  EXPECT_EQ(records[1].text, "before 42 after");
  EXPECT_EQ(records[2].text, "inner");
  EXPECT_EQ(records[3].text, "before 42 after");
}

// Lines of two objects on the same thread are passed on separately
TEST_F(LoggingTest, ObjectsStreamSeparateLines)
{
  const Logger a;
  const Logger b;
  a.part("A");
  EXPECT_TRUE(m_sink->records().empty());
  b.part("B");  // ends the line of a
  b.end();
  const auto records = m_sink->records();
  ASSERT_EQ(records.size(), 2U);
  EXPECT_EQ(records[0].text, "A");
  EXPECT_EQ(records[1].text, "B");
}

// Threads that stream through one shared object each build their own line
TEST_F(LoggingTest, ThreadsSharingAnObject)
{
  const Logger log;
  constexpr int n_threads = 8;
  constexpr int n_lines = 200;
  {
    std::vector<std::jthread> threads;
    for (int t = 0; t < n_threads; ++t)
      threads.emplace_back(
          [&log]
          {
            for (int i = 0; i < n_lines; ++i)
              log.streamed(5);
          });
  }
  const auto records = m_sink->records();
  ASSERT_EQ(records.size(), std::size_t {n_threads * n_lines});
  for (const auto& record : records)
    EXPECT_EQ(record.text, "items: 0 1 2 3 4");
}

// Sinks can be replaced while other threads log, and no message is lost
TEST(LoggingSinkTest, ReplaceWhileLogging)
{
  const auto first = std::make_shared<CaptureSink>();
  const auto second = std::make_shared<CaptureSink>();
  FastCaloSim::set_log_sink(first);
  const Logger log;
  constexpr int n_threads = 4;
  constexpr int n_messages = 2000;
  {
    std::vector<std::jthread> threads;
    for (int t = 0; t < n_threads; ++t)
      threads.emplace_back(
          [&log]
          {
            for (int i = 0; i < n_messages; ++i)
              log.info(i);
          });
    for (int i = 0; i < 1000; ++i)
      FastCaloSim::set_log_sink(i % 2 == 0 ? second : first);
  }
  FastCaloSim::set_log_sink(nullptr);
  EXPECT_EQ(first->records().size() + second->records().size(),
            std::size_t {n_threads * n_messages});
}

// A failing sink neither throws at the call site nor loses the message
TEST(LoggingSinkTest, ThrowingSinkFallsBackToDefault)
{
  FastCaloSim::set_log_sink(std::make_shared<ThrowingSink>());
  const Logger log;
  testing::internal::CaptureStdout();
  EXPECT_NO_THROW(log.info(5));
  const std::string output = testing::internal::GetCapturedStdout();
  FastCaloSim::set_log_sink(nullptr);
  EXPECT_TRUE(output.ends_with(" (INFO) info 5\n"));
}

// Destructors of static objects can still stream lines after the thread-local
// objects of the main thread are destroyed
TEST(LoggingDeathTest, StreamFromStaticDestructor)
{
  EXPECT_EXIT(
      {
        FastCaloSim::set_log_sink(std::make_shared<CerrSink>());
        const Logger log;
        log.streamed(1);
        static const ExitLogger exit_logger(true);
        std::exit(0);
      },
      ::testing::ExitedWithCode(0),
      "streamed at exit");
}

// Also when nothing ends the line
TEST(LoggingDeathTest, UnfinishedLineFromStaticDestructor)
{
  EXPECT_EXIT(
      {
        FastCaloSim::set_log_sink(std::make_shared<CerrSink>());
        const Logger log;
        log.streamed(1);
        static const ExitLogger exit_logger(false);
        std::exit(0);
      },
      ::testing::ExitedWithCode(0),
      "streamed at exit");
}

TEST(LoggingDefaultSinkTest, WritesDecoratedLine)
{
  const Logger log;
  testing::internal::CaptureStdout();
  log.info(3);
  const std::string output = testing::internal::GetCapturedStdout();
  EXPECT_TRUE(output.starts_with("   /LoggingTests.cxx:"));
  EXPECT_TRUE(output.ends_with(" (INFO) info 3\n"));
  EXPECT_EQ(std::ranges::count(output, '\n'), 1);
}
