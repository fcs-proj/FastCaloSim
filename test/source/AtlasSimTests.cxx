#include "AtlasSimTests.h"

#include <gtest/gtest.h>

#include "FastSimModel.hh"
#include "G4GlobalFastSimulationManager.hh"
#include "TestConfig/AtlasSimTestsConfig.h"
#include "TestHelpers/Event.h"
#include "TestHelpers/IOManager.h"
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
  model->setExtrapolationTool(extrapolator);

  // Activate the simulation for this test
  model->doSimulation(true);

  // Set the param file path
  std::string paramsObject {"SelPDGID"};
  TFCSParametrizationBase* param = static_cast<TFCSParametrizationBase*>(
      BasicSimTests::param_file->Get(paramsObject.c_str()));

  // Set up the geometry
  param->set_geometry(AtlasGeoTests::geo);
  // Set logging level
  param->setLevel(MSG::Level::VERBOSE);

  // Set the param in the model
  model->setParametrization(param);

  // Start a run
  G4RunTests::run_manager->BeamOn(1);

  // Retrieve the vector of simulation states
  TestHelpers::SimStateContainer states = model->GetSimulationStates();

  // Serialize simulation states to json file
  states.serialize(AtlasSimTestConfig::sim_output_path());

  // Plot the cells containing energy and colour according to the energy
  std::string sim_plot_exec = AtlasSimTestConfig::py_sim_exec(evt.get_label());

  ASSERT_TRUE(system(sim_plot_exec.c_str()) == 0);

  // Total simulated energy
  ASSERT_NE(states.at(0).E(), 32070.4);
  // Simulated energy in sampling 1
  ASSERT_NE(states.at(0).E(1), 11161.9);
}

INSTANTIATE_TEST_SUITE_P(AtlasSimulation,
                         AtlasSimTests,
                         ::testing::ValuesIn(AtlasSimTestConfig::EVENTS),
                         testing::PrintToStringParamName());
