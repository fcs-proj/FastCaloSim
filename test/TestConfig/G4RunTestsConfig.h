#pragma once

#include <gtest/gtest.h>

class G4RunTestsConfig : public ::testing::Environment
{
public:
  // The location of the simplified geometry file
  inline static const std::string GEO_FILE_PATH =
      std::string(TEST_BASE_DIR) + "/data/geometry/ATLAS/simplified_geo.gdml";
  // Boolean flag to enable verbose tracking
  inline static const bool ENABLE_VERBOSE_TRACKING = false;
  // Boolean flag to enable usage of the ATLAS magnetic field
  inline static const bool USE_ATLAS_FIELD = true;
};
