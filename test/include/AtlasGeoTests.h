#pragma once

#include <CLHEP/Random/RanluxEngine.h>
#include <gtest/gtest.h>

#include "FastCaloSim/Geometry/CaloCellLookup.h"
#include "ROOT/RDataFrame.hxx"
#include "TestConfig/AtlasGeoTestsConfig.h"
#include "TestHelpers/Event.h"

// Test suite for GeoTests
class AtlasGeoTests
    : public virtual ::testing::TestWithParam<TestHelpers::Event>
{
protected:
  static CaloCellLookup* lookup;

  // Sets up the test suite
  // Called before the first test in this test suite
  static void SetUpTestSuite()
  {
    ROOT::RDataFrame df =
        ROOT::RDataFrame(AtlasGeoTestsConfig::ATLAS_CALO_CELL_TREE_NAME,
                         AtlasGeoTestsConfig::ATLAS_CALO_CELL_PATH);
    lookup = new CaloCellLookup(df);
  }

  // Tears down the test suite
  // Called after the last test in this test suite
  static void TearDownTestSuite() { delete lookup; }
};

// Initialize the static member
CaloCellLookup* AtlasGeoTests::lookup = nullptr;
