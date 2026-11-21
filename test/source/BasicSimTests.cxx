// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#include "BasicSimTests.h"

#include <CLHEP/Random/RanluxEngine.h>
#include <gtest/gtest.h>

#include "FastCaloSim/Core/TFCSExtrapolationState.h"
#include "FastCaloSim/Core/TFCSParametrizationBase.h"
#include "FastCaloSim/Core/TFCSSimulationState.h"
#include "FastCaloSim/Core/TFCSTruthState.h"

TEST_F(BasicSimTests, ReadParamFile)
{
  Long64_t const bytes_read = param_files["barrel"]->GetBytesRead();
  EXPECT_EQ(bytes_read, 10097);
}

TEST_F(BasicSimTests, CheckParamObject)
{
  const std::string paramsObject {"SelPDGID"};
  TFCSParametrizationBase const* param = static_cast<TFCSParametrizationBase*>(
      param_files["barrel"]->Get(paramsObject.c_str()));
  EXPECT_NE(param, nullptr);
}

TEST_F(BasicSimTests, DoDummySimulation)
{
  std::string const paramsObject {"SelPDGID"};
  TFCSParametrizationBase* param = static_cast<TFCSParametrizationBase*>(
      param_files["barrel"]->Get(paramsObject.c_str()));

  // Set up the geometry
  param->set_geometry(AtlasGeoTests::geo);
  // Set logging level
  param->setLevel(FastCaloSim::MSG::Level::VERBOSE);

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
  for (std::size_t i = 0; i < 24; ++i) {
    // r and z increase by 10 mm per layer
    const double offset = 10.0 * static_cast<double>(i);
    extrapol_state.set_eta(i, FastCaloSim::SubPos::ENT, truth_state.Eta());
    extrapol_state.set_eta(i, FastCaloSim::SubPos::EXT, truth_state.Eta());
    extrapol_state.set_eta(i, FastCaloSim::SubPos::MID, truth_state.Eta());
    extrapol_state.set_phi(i, FastCaloSim::SubPos::ENT, truth_state.Phi());
    extrapol_state.set_phi(i, FastCaloSim::SubPos::EXT, truth_state.Phi());
    extrapol_state.set_phi(i, FastCaloSim::SubPos::MID, truth_state.Phi());
    extrapol_state.set_r(i, FastCaloSim::SubPos::ENT, 1500 + offset);
    extrapol_state.set_r(i, FastCaloSim::SubPos::EXT, 1510 + offset);
    extrapol_state.set_r(i, FastCaloSim::SubPos::MID, 1505 + offset);
    extrapol_state.set_z(i, FastCaloSim::SubPos::ENT, 3500 + offset);
    extrapol_state.set_z(i, FastCaloSim::SubPos::EXT, 3510 + offset);
    extrapol_state.set_z(i, FastCaloSim::SubPos::MID, 3505 + offset);
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
