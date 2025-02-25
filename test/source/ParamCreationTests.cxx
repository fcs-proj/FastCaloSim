#include <cstdlib>  // For std::system

#include <gtest/gtest.h>

#include "TestHelpers/IOManager.h"
#include "TestHelpers/LogComparer.h"

TEST(ParamCreationTests, Barrel)
{
  // Create the output directory for the test
  const std::string output_dir =
      TestHelpers::IOManager::create_test_output_dir();

  // Get the path to the reference files
  const std::string ref_dir = TestHelpers::IOManager::ref_dir();

  // Path to output and reference log files
  const std::string log_file = output_dir + "out.log";
  const std::string ref_log_file = ref_dir + "out.log";

  // Build the path to the executable
  const std::string EXECUTABLE_PATH =
      std::string(PARAM_BINARY_DIR) + "createParamSlice";

  // Run the executable
  std::string command = EXECUTABLE_PATH + " > " + log_file + " 2>&1";
  int returnCode = std::system(command.c_str());
  EXPECT_EQ(returnCode, 0) << "Executable failed to run. Check log: "
                           << log_file;

  // Compare the log files
  LogComparer comparer(ref_log_file, log_file);

  ASSERT_TRUE(comparer.compareLogs())
      << "Log files differ! See the filtered differences above.";
}
