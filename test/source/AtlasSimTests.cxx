#include "AtlasSimTests.h"

#include <gtest/gtest.h>

#include "FastSimModel.hh"
#include "G4GlobalFastSimulationManager.hh"
#include "TestConfig/AtlasSimTestsConfig.h"
#include "TestHelpers/Event.h"
#include "TestHelpers/IOManager.h"
#include "TestHelpers/JsonComparer.h"
#include "TestHelpers/SimStateContainer.h"

TEST_P(AtlasSimTests, AtlasSimulation)
{
  // Create the output directory of the test
  const std::string output_dir =
      TestHelpers::IOManager::create_test_output_dir();

  // Set the output directory to the configuration
  AtlasSimTestConfig::set_output_dir(output_dir);

  // Retrieve the event to be transported
  const TestHelpers::Event evt = GetParam();

  // Tell the action initialization to set the particles
  G4RunTests::actionInitialization->set_particle_container(
      &evt.get_container());
  G4RunTests::run_manager->SetUserInitialization(
      G4RunTests::actionInitialization);

  // Initialize G4 kernel
  G4RunTests::run_manager->Initialize();

  /* Configure the fast simulation model */
  auto* mgr = G4GlobalFastSimulationManager::GetGlobalFastSimulationManager();
  auto* model =
      dynamic_cast<FastSimModel*>(mgr->GetFastSimulationModel("model"));

  // Reset the state of the fast simulation model
  model->Reset();

  // Set the extrapolation tool
  BasicExtrapolTests::extrapolator.setLevel(
      AtlasSimTestConfig::EXTRAPOLATION_MSG_LEVEL);
  model->setExtrapolationTool(BasicExtrapolTests::extrapolator);

  // Activate the simulation for this test
  model->doSimulation(true);

  // Set the param file path
  std::string paramsObject {"SelPDGID"};
  TFCSParametrizationBase* param = static_cast<TFCSParametrizationBase*>(
      BasicSimTests::param_files[evt.get_param_file_key()]->Get(
          paramsObject.c_str()));

  // Set up the geometry
  param->set_geometry(AtlasDeprecatedGeoTests::geo);
  // Set logging level
  param->setLevel(MSG::Level::VERBOSE);

  // Set the param in the model
  model->setParametrization(param);

  // Start a run
  G4RunTests::run_manager->BeamOn(1);

  // Retrieve the transportation tracks and serialize to json file
  TestHelpers::TrackContainer tracks = model->GetTransportTracks();
  tracks.serialize(AtlasSimTestConfig::transport_output_path());

  // Retrieve the simulation states and serialize to json file
  TestHelpers::SimStateContainer states = model->GetSimulationStates();
  states.serialize(AtlasSimTestConfig::sim_output_path());

  // Plot the cells containing energy and color according to the energy
  std::string sim_plot_exec = AtlasSimTestConfig::py_sim_exec(evt.get_label());

  ASSERT_TRUE(system(sim_plot_exec.c_str()) == 0);

  // Now compare results with references with a tolerance of 0.1 percent
  const double tol = 1e-3;
  JsonComparer json(tol);

  // Compare the transport data to references
  ASSERT_TRUE(json.compare(AtlasSimTestConfig::transport_output_path(),
                           AtlasSimTestConfig::transport_output_ref_path()));
  // Compare the simulation data to references
  ASSERT_TRUE(json.compare(AtlasSimTestConfig::sim_output_path(),
                           AtlasSimTestConfig::sim_output_ref_path()));
}

INSTANTIATE_TEST_SUITE_P(AtlasSimulation,
                         AtlasSimTests,
                         ::testing::ValuesIn(AtlasSimTestConfig::EVENTS),
                         testing::PrintToStringParamName());
