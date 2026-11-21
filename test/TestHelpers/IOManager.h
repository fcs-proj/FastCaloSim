// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#pragma once

#include <filesystem>

#include <gtest/gtest.h>

namespace TestHelpers::IOManager
{
inline auto test_suite_name() -> std::string
{
  return ::testing::UnitTest::GetInstance()
      ->current_test_info()
      ->test_suite_name();
}

inline auto test_name() -> std::string
{
  return ::testing::UnitTest::GetInstance()->current_test_info()->name();
}

inline auto create_test_output_dir() -> std::string
{
  const std::string output_dir =
      TEST_OUTPUT_DIR + test_suite_name() + "_" + test_name() + "/";
  // Create the output directory if it does not exist
  std::filesystem::create_directories(output_dir);

  return output_dir;
}

inline auto ref_dir() -> std::string
{
  return TEST_REFS_DIR + test_suite_name() + "_" + test_name() + "/";
}

}  // namespace TestHelpers::IOManager
