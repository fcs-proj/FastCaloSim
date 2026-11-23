// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#include <cmath>
#include <stdexcept>
#include <vector>

#include "BasicExtrapolTests.h"

#include <gtest/gtest.h>

#include <G4ThreeVector.hh>

#include "FastCaloSim/Core/TFCSExtrapolationState.h"
#include "FastCaloSim/Core/TFCSTruthState.h"
#include "FastCaloSim/Extrapolation/FastCaloSimCaloExtrapolation.h"
#include "G4FieldTrack.hh"

TEST_F(BasicExtrapolTests, BasicExtrapolation)
{
  // Define a dummy extrapolation state
  TFCSExtrapolationState extrapol_state;
  // Define a dummy truth state
  TFCSTruthState truth_state;
  truth_state.SetPtEtaPhiM(65536, 0.225, 1.8, 0);
  truth_state.set_pdgid(22);

  // Vector of positions to be used in the extrapolation
  std::vector<G4ThreeVector> const pos_vec = {
      G4ThreeVector(1429.33, 64.2977, 2229.77),
      G4ThreeVector(1482.48, -13.6767, 2325.32),
      G4ThreeVector(1502.78, -63.3824, 2359.21),
      G4ThreeVector(1546.13, -269.569, 2403.29),
      G4ThreeVector(1567.91, -696.375, 2004.09),
      G4ThreeVector(1731.87, -645.539, 1357.94),
      G4ThreeVector(1733.41, -644.271, 1354.47),
      G4ThreeVector(1890.88, -504.331, 1063.16),
      G4ThreeVector(2297.71, -102.624, 554.826),
      G4ThreeVector(2594.23, 173.105, 246.015),
      G4ThreeVector(2833.65, 362.46, 2.75335e-14),
      G4ThreeVector(3354.55, 761.972, -465.1),
      G4ThreeVector(3680.91, 1025.15, -639.865)};

  // Dummy values for the momentum direction and polarization
  // Note: the momentum direction only has an influence on
  // angle between extrapolated position and momentum (Angle3D and AngleEta)
  G4ThreeVector const mom_dir = G4ThreeVector(0, 1, 1.0);
  G4ThreeVector const pol = G4ThreeVector(0.0, 0.0, 0.0);

  std::vector<G4FieldTrack> calo_steps;
  calo_steps.reserve(pos_vec.size());
  for (const auto& pos : pos_vec) {
    G4FieldTrack const track = G4FieldTrack(pos, 0, mom_dir, 0, 0.0, 0.0, pol);
    calo_steps.emplace_back(track);
  }

  // Perform the extrapolation
  extrapolator.extrapolate(extrapol_state, &truth_state, calo_steps);

  // Test extrapolation to ID-Calo boundary
  ASSERT_NEAR(extrapol_state.IDCaloBoundary_eta(), 1.41753, 1e-1);
  ASSERT_NEAR(extrapol_state.IDCaloBoundary_phi(), 0.0449542, 1e-1);
  ASSERT_NEAR(extrapol_state.IDCaloBoundary_r(), 1148, 1e-1);

  // Helper lambda to test layer values
  auto test_layer =
      [&](std::size_t layer, double eta, double phi, double r, double z) -> void
  {
    ASSERT_NEAR(extrapol_state.eta(layer, FastCaloSim::SubPos::MID), eta, 1e-1);
    ASSERT_NEAR(extrapol_state.phi(layer, FastCaloSim::SubPos::MID), phi, 1e-1);
    ASSERT_NEAR(extrapol_state.r(layer, FastCaloSim::SubPos::MID), r, 1e-1);
    ASSERT_NEAR(extrapol_state.z(layer, FastCaloSim::SubPos::MID), z, 1e-1);
  };

  // Test each layerß
  test_layer(0, 1.22971, 0.018368, 1455.18, 2275.81);
  test_layer(1, 1.22919, -0.0691782, 1515.1, 2368.05);
  test_layer(2, 0.96349, -0.411096, 1729.23, 1936.11);
  test_layer(3, 0.641789, -0.332594, 1872.92, 1286.25);
  test_layer(4, 1.01214, -0.402097, 3072.79, 3669);
  test_layer(5, 0.454442, -0.207122, 8061.18, 3790.74);
  test_layer(6, 0.337639, -0.1168, 11757.8, 4045.74);
  test_layer(7, 0.515638, -0.25734, 7831.9, 4219.78);
  test_layer(8, 0.408492, -0.170821, 10623.3, 4461.25);
  test_layer(9, 0.26116, -0.0606182, 18435.4, 4869.5);
  test_layer(10, 0.261947, -0.0611869, 20473.4, 5424.5);
  test_layer(11, 0.296685, -0.0864763, 19614.3, 5905);
  test_layer(12, 0.161155, 0.0157981, 2450, 396.542);
  test_layer(13, 0.0208621, 0.113869, 2795, 58.3139);
  test_layer(14, 1.23227, -0.0295309, 3630, 5694.34);
  test_layer(15, -0.0895637, 0.190935, 3215, -288.332);
  test_layer(16, -0.151687, 0.248825, 3630, -552.737);
  test_layer(17, -0.0321279, 0.150037, 110041.23, -3536);
  test_layer(18, 0.161155, 0.0157981, 2450, 396.542);
  test_layer(19, -0.00378506, 0.129909, 2870, -10.8632);
  test_layer(20, -0.138533, 0.228975, 3480, -483.639);
  test_layer(21, 2.2997, 0.0449542, 1000.03, 4935.6);
  test_layer(22, 2.2869, 0.0449542, 1107.37, 5394.45);
  test_layer(23, 2.29705, 0.0449542, 1192.37, 5868.95);

  // // Check the angles
  ASSERT_NEAR(extrapol_state.IDCaloBoundary_AngleEta(), -0.30994, 1e-1);
  ASSERT_NEAR(extrapol_state.IDCaloBoundary_Angle3D(), 0.872097, 1e-1);
}

// A particle that enters the calorimeter from the beam pipe, through the FCal
// bore beyond the end of the ID-calo boundary (here a pion of a ttbar event
// in Athena): it reaches the boundary on the endcap of the innermost cylinder,
// which is at the same |z| as the end of the previous one
TEST_F(BasicExtrapolTests, ExtrapolationFromBeamPipe)
{
  constexpr double eta = 3.5396146600687683;
  constexpr double phi = 1.0748289639304982;
  for (const double side : {1.0, -1.0}) {
    TFCSExtrapolationState extrapol_state;
    TFCSTruthState truth_state;
    truth_state.SetPtEtaPhiM(800.9999963160438, side * eta, phi, 139.57);
    truth_state.set_pdgid(211);
    truth_state.set_vertex(
        14.419529402308543, 38.380687484931265, side * 4683.227376579961);

    G4ThreeVector direction;
    direction.setRThetaPhi(1, 2 * std::atan(std::exp(-side * eta)), phi);
    G4ThreeVector const start(truth_state.vertex().X(),
                              truth_state.vertex().Y(),
                              truth_state.vertex().Z());
    std::vector<G4FieldTrack> calo_steps;
    for (const double distance : {0.0, 100.0, 200.0}) {
      calo_steps.emplace_back(start + distance * direction,
                              0,
                              direction,
                              0,
                              0.0,
                              0.0,
                              G4ThreeVector());
    }

    extrapolator.extrapolate(extrapol_state, &truth_state, calo_steps);

    ASSERT_NE(extrapol_state.IDCaloBoundary_eta(), -999.);
    EXPECT_NEAR(extrapol_state.IDCaloBoundary_z(), side * 4587, 1e-3);
    // Straight back along the track to |z| = 4587 mm
    EXPECT_NEAR(extrapol_state.IDCaloBoundary_r(), 35.4696, 1e-3);
    EXPECT_NEAR(extrapol_state.IDCaloBoundary_phi(), 1.23288, 1e-4);
  }
}

// A track that starts just beyond the endcap of the outer boundary cylinder
// leaves the ID where it crosses that endcap, not at the point of the inner
// cylinder closest to the track, which it never reaches but is closer to the
// origin
TEST_F(BasicExtrapolTests, ExtrapolationJustBeyondEndcap)
{
  for (const double side : {1.0, -1.0}) {
    TFCSExtrapolationState extrapol_state;
    TFCSTruthState truth_state;
    truth_state.SetPtEtaPhiM(10000, side * 2.6, 0.38, 0);
    truth_state.set_pdgid(22);

    const G4ThreeVector direction = G4ThreeVector(1, 0, side).unit();
    const G4ThreeVector start(500, 200, side * 3555);
    std::vector<G4FieldTrack> calo_steps;
    for (const double distance : {0.0, 100.0, 200.0}) {
      calo_steps.emplace_back(start + distance * direction,
                              0,
                              direction,
                              0,
                              0.0,
                              0.0,
                              G4ThreeVector());
    }

    extrapolator.extrapolate(extrapol_state, &truth_state, calo_steps);

    ASSERT_NE(extrapol_state.IDCaloBoundary_eta(), -999.);
    // Straight back along the track to |z| = 3550 mm, at x = 495 mm
    EXPECT_NEAR(extrapol_state.IDCaloBoundary_z(), side * 3550, 1e-3);
    EXPECT_NEAR(extrapol_state.IDCaloBoundary_x(), 495, 1e-3);
    EXPECT_NEAR(extrapol_state.IDCaloBoundary_y(), 200, 1e-3);
  }
}

// The extrapolation has no detector of its own: the host must set it up
TEST_F(BasicExtrapolTests, ExtrapolationNeedsSettings)
{
  TFCSExtrapolationState extrapol_state;
  TFCSTruthState truth_state;
  truth_state.SetPtEtaPhiM(65536, 0.225, 1.8, 0);
  const std::vector<G4FieldTrack> calo_steps;

  FastCaloSimCaloExtrapolation without_settings;
  without_settings.set_geometry(AtlasGeoTests::geo);
  EXPECT_THROW(
      without_settings.extrapolate(extrapol_state, &truth_state, calo_steps),
      std::logic_error);

  FastCaloSimCaloExtrapolation without_geometry;
  without_geometry.set_settings(atlasSettings());
  EXPECT_THROW(
      without_geometry.extrapolate(extrapol_state, &truth_state, calo_steps),
      std::logic_error);
}

TEST_F(BasicExtrapolTests, InconsistentSettings)
{
  FastCaloSimCaloExtrapolation tool;

  auto settings = atlasSettings();
  settings.id_calo_boundary.clear();
  EXPECT_THROW(tool.set_settings(settings), std::invalid_argument);

  // Each cylinder must be closer to the beam axis than the previous one
  settings = atlasSettings();
  settings.id_calo_boundary[1].radius = 1148.0;
  EXPECT_THROW(tool.set_settings(settings), std::invalid_argument);

  // ... and reach at least as far in |z|
  settings = atlasSettings();
  settings.id_calo_boundary[1].half_length = 3000.0;
  EXPECT_THROW(tool.set_settings(settings), std::invalid_argument);

  // Only barrel layers can end at another layer, and only at one the geometry
  // has, checked whichever of settings and geometry is set last
  tool.set_geometry(AtlasGeoTests::geo);
  settings = atlasSettings();
  settings.barrel_layer_ends_at_layer[5] = 6;  // EME1 is not a barrel layer
  EXPECT_THROW(tool.set_settings(settings), std::invalid_argument);
  settings = atlasSettings();
  settings.barrel_layer_ends_at_layer[0] = 1000;
  EXPECT_THROW(tool.set_settings(settings), std::invalid_argument);

  FastCaloSimCaloExtrapolation geometry_last;
  geometry_last.set_settings(settings);
  EXPECT_THROW(geometry_last.set_geometry(AtlasGeoTests::geo),
               std::invalid_argument);

  EXPECT_NO_THROW(tool.set_settings(atlasSettings()));
}
