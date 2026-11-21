// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

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
                 "TFCSParam_pid22_pid211_Mom16384_131072_eta_315_325.root"},
          // Slice of the ATLAS MC23 parametrization
          // (TFCSparam_dev_Hybrid_Ha_v5_all_baryons_0_500.root), all energies
          // of one barrel eta bin, objects as written by the generation. It
          // covers the GAN, the ML extrapolation weights and the pion shower
          // shapes.
          {"barrel_ml",
           std::string(TEST_BASE_DIR)
               + "/data/param/barrel_ml/"
                 "TFCSParam_pid22_pid211_Mom64_4194304_eta_20_25.root"}};

  inline static auto param_file_path(const std::string& region) -> std::string
  {
    return PARAM_FILE_PATHS.at(region);
  }
};
