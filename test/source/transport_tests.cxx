#include <gtest/gtest.h>

#include "ActionInitialization.hh"
#include "DetectorConstruction.hh"
#include "FTFP_BERT.hh"
#include "FastSimModel.hh"
#include "G4FastSimulationPhysics.hh"
#include "G4GDMLParser.hh"
#include "G4GlobalFastSimulationManager.hh"
#include "G4RunManagerFactory.hh"
#include "TestHelpers/Event.h"
#include "TestHelpers/IOManager.h"
#include "TestHelpers/ParticleContainer.h"
#include "TestHelpers/ParticleSampler.h"
#include "TestHelpers/ParticleTypes.h"

class TestEnvironment : public ::testing::Environment
{
public:
  // The location of the python script to plot the transport
  inline static const std::string PYTHON_SCRIPT =
      std::string(PYTHON_SCRIPT_DIR) + "plot_transport.py";
  // The location of calo layer data used for plotting the transport
  inline static const std::string CALO_LAYER_DATA =
      std::string(PYTHON_SCRIPT_DIR) + "calo_layers.csv";
  // The location of the simplified geometry file
  inline static const std::string GEO_FILE_PATH =
      std::string(TEST_GEO_DIR) + "simplified_geo.gdml";
  // Boolean flag to enable verbose tracking
  inline static const bool ENABLE_VERBOSE_TRACKING = false;
  // Boolean flag to enable usage of the ATLAS magnetic field
  inline static const bool USE_ATLAS_FIELD = true;
  // Boolean flag to zoom into region with tracks in plotting
  inline static const bool ZOOM_IN = false;

  // The list of events (event = particle container) to be transported
  inline static const std::vector<TestHelpers::Event> EVENT_VECTOR = []()
  {
    // Vector of events to be processed
    std::vector<TestHelpers::Event> events;

    // For each event we sample particles of the same type uniformly in eta
    const float min_eta = -5;
    const float max_eta = 5;
    const float step_size = 0.5;
    const float energy = 50 * CLHEP::MeV;
    std::vector<TestHelpers::ParticleType> ptype = {
        TestHelpers::ParticleTypes::Photon,
        TestHelpers::ParticleTypes::Electron,
        TestHelpers::ParticleTypes::Positron,
        TestHelpers::ParticleTypes::Proton};

    // Particle sampler on the ID-calorimeter boundary
    TestHelpers::IDCaloBoundarySampler sampler;

    for (const auto& type : ptype) {
      TestHelpers::ParticleContainer particles = sampler.uniformEtaSample(
          type.pid, energy, min_eta, max_eta, step_size);

      // energy string rounded to 2 decimal places
      std::stringstream estring;
      estring << std::fixed << std::setprecision(2) << energy;
      std::string evt_label = "$ E=" + estring.str() + R"(\,\text{MeV}\,)"
          + std::string(type.label) + "$";
      // Construct the event with the particle container and event label
      TestHelpers::Event evt(particles, evt_label);

      events.emplace_back(evt);
    }

    return events;
  }();
};

class TransportTests : public ::testing::TestWithParam<TestHelpers::Event>
{
protected:
  static G4RunManager* runManager;

  // Per-test-suite set-up.
  // Called before the first test in this test suite.
  static void SetUpTestSuite()
  {
    // Load the simplified geometry from a GDML file into the volume store
    G4GDMLParser fParser;
    fParser.Read(TestEnvironment::GEO_FILE_PATH, false);

    // Construct the default run manager
    runManager =
        G4RunManagerFactory::CreateRunManager(G4RunManagerType::Serial);

    // Construct the detector
    auto* detectorConstruction = new DetectorConstruction();
    detectorConstruction->setUseAtlasField(TestEnvironment::USE_ATLAS_FIELD);
    runManager->SetUserInitialization(detectorConstruction);

    // Set verbosity of the tracking manager
    G4TrackingManager* trackingManager =
        G4EventManager::GetEventManager()->GetTrackingManager();
    trackingManager->SetVerboseLevel(
        static_cast<G4int>(TestEnvironment::ENABLE_VERBOSE_TRACKING));

    // Configure fast simulation
    auto* fastSimulationPhysics = new G4FastSimulationPhysics();
    fastSimulationPhysics->BeVerbose();

    // Activate fast simulation for the following particles
    std::vector<G4String> fast_sim_particles = {
        "e-", "e+", "gamma", "pi+", "pi-", "pi0", "proton", "anti-proton"};
    for (auto& particle : fast_sim_particles) {
      fastSimulationPhysics->ActivateFastSimulation(particle);
    }

    // Create physics list with fast simulation
    auto* physicsList = new FTFP_BERT(0);
    physicsList->RegisterPhysics(fastSimulationPhysics);
    runManager->SetUserInitialization(physicsList);
  }

  // Per-test-suite tear-down.
  // Called after the last test in this test suite.
  static void TearDownTestSuite()
  {
    // Terminate run
    delete runManager;
  }
};

// Define static members to make them accessible in the tests
G4RunManager* TransportTests::runManager = nullptr;

TEST_P(TransportTests, ParticleTransportTest)
{
  // Retrieve the particle to be transported
  const TestHelpers::Event evt = GetParam();
  // Create the output directory of the test
  const std::string output_dir =
      TestHelpers::IOManager::create_test_output_dir();

  // Create action initialization to set up particle to simulate
  auto* actionInitialization = new ActionInitialization();

  // Tell the action initialization to set the particl
  actionInitialization->set_particle_container(&evt.get_container());
  runManager->SetUserInitialization(actionInitialization);

  // Initialize G4 kernel
  runManager->Initialize();

  // Set the path where the transport vector will be saved to in the fast
  // simulation model
  const std::string transport_output_path =
      output_dir + "particle_transport.csv";
  auto* mgr = G4GlobalFastSimulationManager::GetGlobalFastSimulationManager();
  auto* model =
      dynamic_cast<FastSimModel*>(mgr->GetFastSimulationModel("model"));
  model->saveTransport(transport_output_path);

  // Start a run
  runManager->BeamOn(1);

  // Set the path where the particle transport will be plotted to
  std::string plot_output_path = output_dir + "particle_transport.png";

  // Run the python script to plot the transport and check that it runs
  // successfully
  std::string command = "python3 " + TestEnvironment::PYTHON_SCRIPT
      + " --input " + transport_output_path + " --output " + plot_output_path
      + " --calo_layer_csv " + TestEnvironment::CALO_LAYER_DATA + " --label \""
      + evt.get_label() + "\"";

  if (TestEnvironment::ZOOM_IN) {
    command += " --track_zoom";
  }

  ASSERT_TRUE(system(command.c_str()) == 0);
}

INSTANTIATE_TEST_SUITE_P(
    ParticleTransportTests,
    TransportTests,
    ::testing::ValuesIn(TestEnvironment::EVENT_VECTOR),
    [](const testing::TestParamInfo<TransportTests::ParamType>& info)
        -> std::string
    {
      std::ostringstream name;
      name << "Event_" << info.index;
      return name.str();
    });