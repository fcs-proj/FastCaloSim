// Copyright (c) 2024 CERN for the benefit of the FastCaloSim project

#include "AtlasTransportTests.h"

#include <gtest/gtest.h>

#include "FastSimModel.hh"
#include "G4GlobalFastSimulationManager.hh"
#include "TestConfig/AtlasTransportTestsConfig.h"
#include "TestHelpers/Event.h"
#include "TestHelpers/IOManager.h"
#include "TestHelpers/JsonComparer.h"

TEST_P(AtlasTransportTests, AtlasTransport)
{
  // Create the output directory of the test
  const std::string output_dir =
      TestHelpers::IOManager::create_test_output_dir();

  // Set the output directory to the configuration
  AtlasTransportTestConfig::set_output_dir(output_dir);

  // Retrieve the event to be transported
  const TestHelpers::Event evt = GetParam();

  // Tell the action initialization to set the particles
  actionInitialization->set_particle_container(&evt.get_container());
  run_manager->SetUserInitialization(actionInitialization);

  // Initialize G4 kernel
  run_manager->Initialize();

  /* Configure the fast simulation model */
  auto* mgr = G4GlobalFastSimulationManager::GetGlobalFastSimulationManager();
  auto* model =
      dynamic_cast<FastSimModel*>(mgr->GetFastSimulationModel("model"));

  // Reset the state of the fast simulation model
  model->Reset();

  // Deactivate the simulation for this test
  model->doSimulation(false);

  // Set the extrapolation tool
  BasicExtrapolTests::extrapolator.setLevel(
      AtlasTransportTestConfig::EXTRAPOLATION_FCS_MSG_LEVEL);
  model->setExtrapolationTool(BasicExtrapolTests::extrapolator);

  // Start a run
  run_manager->BeamOn(1);

  // Serialize tracks to json file
  TestHelpers::TrackContainer tracks = model->GetTransportTracks();
  tracks.serialize(AtlasTransportTestConfig::transport_output_path());

  // Construct the python command to plot the transport
  std::string transport_exec =
      AtlasTransportTestConfig::py_transport_exec(evt.get_label());

  // Plot the transport
  ASSERT_TRUE(system(transport_exec.c_str()) == 0);

  // Serialize the extrapolation states to json file
  TestHelpers::ExtrapStateContainer extrapolations = model->GetExtrapolations();
  extrapolations.serialize(
      AtlasTransportTestConfig::extrapolation_output_path());

  // Construct the python command to plot the extrapolation
  std::string extrapol_exec =
      AtlasTransportTestConfig::py_extrapol_exec(evt.get_label());

  // Plot the extrapolation
  ASSERT_TRUE(system(extrapol_exec.c_str()) == 0);

  // Now compare results with references with a tolerance of 0.1 percent
  const double tol = 1e-3;
  JsonComparer json(tol);

  // Compare the transport data to references
  ASSERT_TRUE(
      json.compare(AtlasTransportTestConfig::transport_output_path(),
                   AtlasTransportTestConfig::transport_output_ref_path()));
  // Compare the extrapolation data to references
  ASSERT_TRUE(
      json.compare(AtlasTransportTestConfig::extrapolation_output_path(),
                   AtlasTransportTestConfig::extrapolation_output_ref_path()));
}

INSTANTIATE_TEST_SUITE_P(AtlasTransport,
                         AtlasTransportTests,
                         ::testing::ValuesIn(AtlasTransportTestConfig::EVENTS),
                         testing::PrintToStringParamName());
