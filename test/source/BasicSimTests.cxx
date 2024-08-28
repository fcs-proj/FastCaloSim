#include "BasicSimTests.h"

#include <CLHEP/Random/RanluxEngine.h>
#include <gtest/gtest.h>

#include "FastCaloSim/Core/TFCSExtrapolationState.h"
#include "FastCaloSim/Core/TFCSParametrizationBase.h"
#include "FastCaloSim/Core/TFCSSimulationState.h"
#include "FastCaloSim/Core/TFCSTruthState.h"

TEST_F(BasicSimTests, ReadParamFile)
{
  Long64_t bytes_read = param_files["barrel"]->GetBytesRead();
  EXPECT_EQ(bytes_read, 10097);
}

TEST_F(BasicSimTests, CheckParamObject)
{
  const std::string paramsObject {"SelPDGID"};
  TFCSParametrizationBase* param = static_cast<TFCSParametrizationBase*>(
      param_files["barrel"]->Get(paramsObject.c_str()));
  EXPECT_NE(param, nullptr);
}

TEST_F(BasicSimTests, DoDummySimulation)
{
  std::string paramsObject {"SelPDGID"};
  TFCSParametrizationBase* param = static_cast<TFCSParametrizationBase*>(
      param_files["barrel"]->Get(paramsObject.c_str()));

  // Set up the geometry
  param->set_geometry(AtlasGeoTests::geo);
  // Set logging level
  param->setLevel(MSG::Level::VERBOSE);

  // Set up a random engine
  std::unique_ptr<CLHEP::RanluxEngine> rnd_engine =
      std::make_unique<CLHEP::RanluxEngine>();
  rnd_engine->setSeed(42);

  // Set up simulation state and assign random engine
  TFCSSimulationState simul_state;
  simul_state.setRandomEngine(rnd_engine.get());

  // Set up particle to simulate
  TFCSTruthState truth_state;
  truth_state.SetPtEtaPhiM(65536, 0.225, 1.8, 0);
  truth_state.set_pdgid(22);

  // Set up a dummy extrapolation state
  TFCSExtrapolationState extrapol_state;
  extrapol_state.set_IDCaloBoundary_eta(truth_state.Eta());
  for (int i = 0; i < 24; ++i) {
    extrapol_state.set_eta(
        i, TFCSExtrapolationState::SUBPOS_ENT, truth_state.Eta());
    extrapol_state.set_eta(
        i, TFCSExtrapolationState::SUBPOS_EXT, truth_state.Eta());
    extrapol_state.set_eta(
        i, TFCSExtrapolationState::SUBPOS_MID, truth_state.Eta());
    extrapol_state.set_phi(
        i, TFCSExtrapolationState::SUBPOS_ENT, truth_state.Phi());
    extrapol_state.set_phi(
        i, TFCSExtrapolationState::SUBPOS_EXT, truth_state.Phi());
    extrapol_state.set_phi(
        i, TFCSExtrapolationState::SUBPOS_MID, truth_state.Phi());
    extrapol_state.set_r(i, TFCSExtrapolationState::SUBPOS_ENT, 1500 + i * 10);
    extrapol_state.set_r(i, TFCSExtrapolationState::SUBPOS_EXT, 1510 + i * 10);
    extrapol_state.set_r(i, TFCSExtrapolationState::SUBPOS_MID, 1505 + i * 10);
    extrapol_state.set_z(i, TFCSExtrapolationState::SUBPOS_ENT, 3500 + i * 10);
    extrapol_state.set_z(i, TFCSExtrapolationState::SUBPOS_EXT, 3510 + i * 10);
    extrapol_state.set_z(i, TFCSExtrapolationState::SUBPOS_MID, 3505 + i * 10);
  }

  // Perform the simulation call
  param->simulate(simul_state, &truth_state, &extrapol_state);

  // Check that simulation was successful
  EXPECT_NEAR(simul_state.E(), 66943.7, 1e-1);
  EXPECT_NEAR(simul_state.E(0), 1096.21, 1e-1);
  EXPECT_NEAR(simul_state.E(1), 19513.9, 1e-1);
  EXPECT_NEAR(simul_state.E(2), 46197.6, 1e-1);
  EXPECT_NEAR(simul_state.E(3), 85.8894, 1e-1);
  EXPECT_NEAR(simul_state.E(12), 50.1765, 1e-1);
}
