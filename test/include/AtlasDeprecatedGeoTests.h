#pragma once

#include <CLHEP/Random/RanluxEngine.h>
#include <gtest/gtest.h>

#include "FastCaloSim/Geometry/CaloGeometryFromFile.h"
#include "TestConfig/AtlasDeprecatedGeoTestsConfig.h"
#include "TestHelpers/Event.h"

class AtlasDeprecatedGeoTests
    : public virtual ::testing::TestWithParam<TestHelpers::Event>
{
protected:
  static CaloGeometryFromFile* geo;

  // Set up the test suite
  // Called before the first test in this test suite.
  static void SetUpTestSuite()
  {
    geo = new CaloGeometryFromFile();

    // Load geometry
    geo->LoadGeometryFromFile(AtlasDeprecatedGeoTestsConfig::GEO_PATH,
                              AtlasDeprecatedGeoTestsConfig::GEO_TAG,
                              AtlasDeprecatedGeoTestsConfig::HASH_FILE_PATH);

    // Load FCal geometry from files
    geo->LoadFCalGeometryFromFiles(
        {AtlasDeprecatedGeoTestsConfig::FCAL1_GEO_PATH,
         AtlasDeprecatedGeoTestsConfig::FCAL2_GEO_PATH,
         AtlasDeprecatedGeoTestsConfig::FCAL3_GEO_PATH});
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
CaloGeometryFromFile* AtlasDeprecatedGeoTests::geo = nullptr;
