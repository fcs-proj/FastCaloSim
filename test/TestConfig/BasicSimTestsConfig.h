#pragma once

#include <string>
#include <unordered_map>

#include <gtest/gtest.h>

class BasicSimTestsConfig : public ::testing::Environment
{
public:
  // Mapping of regions to their respective parameter file paths
  inline static const std::unordered_map<std::string, std::string>
      PARAM_FILE_PATHS = {
          {"barrel",
           std::string(TEST_BASE_DIR)
               + "/data/param/barrel/"
                 "TFCSParam_pid22_pid211_Mom16384_131072_eta_15_30.root"},
          {"endcap",
           std::string(TEST_BASE_DIR)
               + "/data/param/endcap/"
                 "TFCSParam_pid22_pid211_Mom16384_131072_eta_195_205.root"},
          {"fcal",
           std::string(TEST_BASE_DIR)
               + "/data/param/fcal/"
                 "TFCSParam_pid22_pid211_Mom16384_131072_eta_345_355.root"},
          {"barrel_endcap_transition",
           std::string(TEST_BASE_DIR)
               + "/data/param/barrel_endcap_transition/"
                 "TFCSParam_pid22_pid211_Mom16384_131072_eta_140_150.root"},
          {"endcap_fcal_transition",
           std::string(TEST_BASE_DIR)
               + "/data/param/endcap_fcal_transition/"
                 "TFCSParam_pid22_pid211_Mom16384_131072_eta_315_325.root"}};

  inline static auto param_file_path(const std::string& region) -> std::string
  {
    return PARAM_FILE_PATHS.at(region);
  }
};
