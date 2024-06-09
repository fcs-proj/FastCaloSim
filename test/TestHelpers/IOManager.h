#pragma once

#include <gtest/gtest.h>

namespace TestHelpers::IOManager
{
inline auto create_test_output_dir() -> std::string
{
  const std::string test_suite_name = ::testing::UnitTest::GetInstance()
                                          ->current_test_info()
                                          ->test_suite_name();
  const std::string test_name =
      ::testing::UnitTest::GetInstance()->current_test_info()->name();

  const std::string output_dir =
      TEST_OUTPUT_DIR + test_suite_name + "_" + test_name + "/";
  // Create the output directory if it does not exist
  const std::string command = "mkdir -p " + output_dir;
  system(command.c_str());

  return output_dir;
}
}  // namespace TestHelpers::IOManager
