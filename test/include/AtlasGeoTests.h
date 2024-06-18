#pragma once

#include <CLHEP/Random/RanluxEngine.h>
#include <gtest/gtest.h>

#include "FastCaloSim/Geometry/CaloGeometryFromFile.h"
#include "TestConfig/AtlasGeoTestsConfig.h"
#include "TestHelpers/Event.h"

class AtlasGeoTests
    : public virtual ::testing::TestWithParam<TestHelpers::Event>
{
protected:
  static CaloGeometryFromFile* geo;

  // Set up the test suite
  // Called before the first test in this test suite.
  // Note: This does not work if incurred via ctest
  static void SetUpTestSuite()
  {
    geo = new CaloGeometryFromFile();

    // Load geometry
    geo->LoadGeometryFromFile(AtlasGeoTestsConfig::GEO_PATH,
                              AtlasGeoTestsConfig::GEO_TAG,
                              AtlasGeoTestsConfig::HASH_FILE_PATH);

    // Load FCal geometry from files
    geo->LoadFCalGeometryFromFiles({AtlasGeoTestsConfig::FCAL1_GEO_PATH,
                                    AtlasGeoTestsConfig::FCAL2_GEO_PATH,
                                    AtlasGeoTestsConfig::FCAL3_GEO_PATH});
  }

  // Tear down the test suite
  // Called after the last test in this test suite.
  static void TearDownTestSuite()
  {
    // Do not delete geo here, as it will lead to crashes when using with
    // AtlasSimTests
  }
};

// Initialize the static members
CaloGeometryFromFile* AtlasGeoTests::geo = nullptr;
