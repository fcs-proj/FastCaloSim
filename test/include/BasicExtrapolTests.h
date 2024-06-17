#pragma once

#include <gtest/gtest.h>

#include "AtlasGeoTests.h"
#include "FastCaloSim/Extrapolation/FastCaloSimCaloExtrapolation.h"

class BasicExtrapolTests : public virtual AtlasGeoTests
{
protected:
  static FastCaloSimCaloExtrapolation extrapolator;
  // Per-test-suite set-up.
  // Called before the first test in this test suite.
  // Note: this does not work if incurred via ctest
  static void SetUpTestSuite()
  {
    // Call the base setup that sets up the geometry
    AtlasGeoTests::SetUpTestSuite();
    // Set the geometry of the extrapolator
    extrapolator.set_geometry(geo);
    // Set the verbosity level
    extrapolator.setLevel(MSG::Level::INFO);
  }

  // Per-test-suite tear-down.
  // Called after the last test in this test suite.
  // Note: this does not work if incurred via ctest
  static void TearDownTestSuite()
  {
    // Call the base tear down that deletes the geometry
    AtlasGeoTests::TearDownTestSuite();
  }
};

// Initialize the static members
FastCaloSimCaloExtrapolation BasicExtrapolTests::extrapolator;