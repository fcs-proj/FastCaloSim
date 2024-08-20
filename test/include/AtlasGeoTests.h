#pragma once

#include <CLHEP/Random/RanluxEngine.h>
#include <gtest/gtest.h>

#include "FastCaloSim/Geometry/CaloGeo.h"
#include "ROOT/RDataFrame.hxx"
#include "TestConfig/AtlasGeoTestsConfig.h"
#include "TestHelpers/Event.h"

// Test suite for GeoTests
class AtlasGeoTests
    : public virtual ::testing::TestWithParam<TestHelpers::Event>
{
protected:
  static CaloGeo* geo;

  // Sets up the test suite
  // Called before the first test in this test suite
  static void SetUpTestSuite()
  {
    ROOT::RDataFrame df =
        ROOT::RDataFrame(AtlasGeoTestsConfig::ATLAS_CALO_CELL_TREE_NAME,
                         AtlasGeoTestsConfig::ATLAS_CALO_CELL_PATH);
    geo = new CaloGeo(df);
  }

  // Tears down the test suite
  // Called after the last test in this test suite
  static void TearDownTestSuite() {}  // deleting geo here causes issues
};

// Initialize the static member
CaloGeo* AtlasGeoTests::geo = nullptr;
