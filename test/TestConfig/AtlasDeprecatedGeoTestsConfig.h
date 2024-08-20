#pragma once

#include <gtest/gtest.h>

class AtlasDeprecatedGeoTestsConfig : public ::testing::Environment
{
public:
  // The location of the ATLAS geometry files
  inline static const std::string GEO_TAG = "ATLAS-R2-2016-01-00-01";
  inline static const std::string GEO_PATH =
      std::string(TEST_BASE_DIR) + "/data/geometry/" + GEO_TAG + ".root";
  inline static const std::string HASH_FILE_PATH = std::string(TEST_BASE_DIR)
      + "/data/geometry/cellId_vs_cellHashId_map.txt";
  inline static const std::string FCAL1_GEO_PATH = std::string(TEST_BASE_DIR)
      + "/data/geometry/FCal1-electrodes.sorted.HV.09Nov2007.dat";
  inline static const std::string FCAL2_GEO_PATH = std::string(TEST_BASE_DIR)
      + "/data/geometry/FCal2-electrodes.sorted.HV.April2011.dat";
  inline static const std::string FCAL3_GEO_PATH = std::string(TEST_BASE_DIR)
      + "/data/geometry/FCal3-electrodes.sorted.HV.09Nov2007.dat";
};
