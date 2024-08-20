#pragma once

#include <gtest/gtest.h>

class BasicSimTestsConfig : public ::testing::Environment
{
public:
  // The location of the param file to be used
  inline static const std::string PARAM_FILE_PATH = std::string(TEST_BASE_DIR)
      + "/data/param/TFCSParam_pid22_pid211_Mom16384_131072_eta_15_30.root";
};
