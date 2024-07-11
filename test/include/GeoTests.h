#pragma once

#include <CLHEP/Random/RanluxEngine.h>
#include <gtest/gtest.h>

#include "ROOT/RDataFrame.hxx"
#include "TestConfig/GeoTestsConfig.h"
#include "TestHelpers/Event.h"

// Test suite for GeoTests
class GeoTests : public virtual ::testing::TestWithParam<TestHelpers::Event>
{
protected:
  static ROOT::RDataFrame geo;

  // Sets up the test suite
  // Called before the first test in this test suite
  static void SetUpTestSuite()
  {
    geo = ROOT::RDataFrame(GeoTestsConfig::ATLAS_CALO_CELL_TREE_NAME,
                           GeoTestsConfig::ATLAS_CALO_CELL_PATH);
  }

  // Tears down the test suite
  // Called after the last test in this test suite
  static void TearDownTestSuite() {}
};

// Initialize the static member
ROOT::RDataFrame GeoTests::geo(GeoTestsConfig::ATLAS_CALO_CELL_TREE_NAME,
                               GeoTestsConfig::ATLAS_CALO_CELL_PATH);
