// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#pragma once

#include <gtest/gtest.h>

#include "AtlasGeoTests.h"
#include "FastCaloSim/Extrapolation/FastCaloSimCaloExtrapolation.h"

class BasicExtrapolTests : public virtual AtlasGeoTests
{
protected:
  static inline FastCaloSimCaloExtrapolation extrapolator;

  /// @brief The extrapolation settings of ATLAS
  static auto atlasSettings() -> FastCaloSimCaloExtrapolation::Settings
  {
    return {
        // The ID-calo boundary (mm)
        .id_calo_boundary = {{.radius = 1148.0, .half_length = 3550.0},
                             {.radius = 120.0, .half_length = 4587.0},
                             {.radius = 41.0, .half_length = 4587.0}},
        // The presampler and the EM barrel layers (PreSamplerB, EMB1-3) end at
        // the front of the EM endcap (EME1); the Tile layers at their own ends
        .barrel_layer_ends_at_layer = {{0, 5}, {1, 5}, {2, 5}, {3, 5}},
        .max_abs_eta = 6,
        .min_warning_pt = 500,
    };
  }
  // Per-test-suite set-up.
  // Called before the first test in this test suite.
  static void SetUpTestSuite()
  {
    // Call the base setup that sets up the geometry
    AtlasGeoTests::SetUpTestSuite();
    // Set the ATLAS settings and the geometry of the extrapolator
    extrapolator.set_settings(atlasSettings());
    extrapolator.set_geometry(AtlasGeoTests::geo);
    // Set the verbosity level
    extrapolator.setLevel(FastCaloSim::MSG::Level::INFO);
  }

  // Per-test-suite tear-down.
  // Called after the last test in this test suite.
  static void TearDownTestSuite()
  {
    // Call the base tear down that deletes the geometry
    AtlasGeoTests::TearDownTestSuite();
  }
};
