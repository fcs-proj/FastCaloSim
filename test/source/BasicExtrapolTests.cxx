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
  std::vector<G4ThreeVector> pos_vec = {
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
  G4ThreeVector mom_dir = G4ThreeVector(0, 1, 1.0);
  G4ThreeVector pol = G4ThreeVector(0.0, 0.0, 0.0);

  std::vector<G4FieldTrack> calo_steps;
  calo_steps.reserve(pos_vec.size());
  for (const auto& pos : pos_vec) {
    G4FieldTrack track = G4FieldTrack(pos, 0, mom_dir, 0, 0.0, 0.0, pol);
    calo_steps.emplace_back(track);
  }

  // Perform the extrapolation
  extrapolator.extrapolate(extrapol_state, &truth_state, calo_steps);

  // Test extrapolation to ID-Calo boundary
  ASSERT_NEAR(extrapol_state.IDCaloBoundary_eta(), 1.41753, 1e-1);
  ASSERT_NEAR(extrapol_state.IDCaloBoundary_phi(), 0.0449542, 1e-1);
  ASSERT_NEAR(extrapol_state.IDCaloBoundary_r(), 1148, 1e-1);

  // Helper lambda to test layer values
  auto test_layer = [&](int layer, double eta, double phi, double r, double z)
  {
    ASSERT_NEAR(
        extrapol_state.eta(layer, TFCSExtrapolationState::SUBPOS::SUBPOS_MID),
        eta,
        1e-1);
    ASSERT_NEAR(
        extrapol_state.phi(layer, TFCSExtrapolationState::SUBPOS::SUBPOS_MID),
        phi,
        1e-1);
    ASSERT_NEAR(
        extrapol_state.r(layer, TFCSExtrapolationState::SUBPOS::SUBPOS_MID),
        r,
        1e-1);
    ASSERT_NEAR(
        extrapol_state.z(layer, TFCSExtrapolationState::SUBPOS::SUBPOS_MID),
        z,
        1e-1);
  };

  // Test each layer
  test_layer(0, 1.22986, 0.0168185, 1456.67, 2278.54);
  test_layer(1, 1.22876, -0.0727093, 1516.63, 2369.21);
  test_layer(2, 0.959902, -0.410334, 1730.76, 1928.51);
  test_layer(3, 0.639429, -0.331143, 1874.45, 1281.94);
  test_layer(4, 1.01214, -0.402097, 3072.79, 3669);
  test_layer(5, 0.454442, -0.207122, 8061.18, 3790.74);
  test_layer(6, 0.337639, -0.1168, 11757.8, 4045.74);
  test_layer(7, 0.515638, -0.25734, 7831.9, 4219.78);
  test_layer(8, 0.408492, -0.170822, 10623.3, 4461.25);
  test_layer(9, 0.26116, -0.0606182, 18435.4, 4869.5);
  test_layer(10, 0.261947, -0.0611871, 20473.4, 5424.5);
  test_layer(11, 0.296685, -0.0864765, 19614.2, 5905);
  test_layer(12, 0.161155, 0.0157981, 2450, 396.542);
  test_layer(13, 0.0208621, 0.113869, 2795, 58.3139);
  test_layer(14, 1.23227, -0.0295309, 3630, 5694.34);
  test_layer(15, -0.0895637, 0.190935, 3215, -288.332);
  test_layer(16, -0.151687, 0.248825, 3630, -552.737);
  test_layer(17, -0.0320941, 0.150013, 110157, -3536);
  test_layer(18, 0.161155, 0.0157981, 2450, 396.542);
  test_layer(19, -0.00378506, 0.129909, 2870, -10.8632);
  test_layer(20, -0.138533, 0.228975, 3480, -483.639);
  test_layer(21, 2.2997, 0.0449542, 1000.03, 4935.6);
  test_layer(22, 2.2869, 0.0449542, 1107.37, 5394.45);
  test_layer(23, 2.29705, 0.0449542, 1192.37, 5868.95);

  // Check the angles
  ASSERT_NEAR(extrapol_state.IDCaloBoundary_AngleEta(), -0.3099, 1e-1);
  ASSERT_NEAR(extrapol_state.IDCaloBoundary_Angle3D(), 0.8721, 1e-1);
}
