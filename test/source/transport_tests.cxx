#include <gtest/gtest.h>

#include "ActionInitialization.hh"
#include "DetectorConstruction.hh"
#include "FTFP_BERT.hh"
#include "FastSimModel.hh"
#include "G4FastSimulationPhysics.hh"
#include "G4GDMLParser.hh"
#include "G4GlobalFastSimulationManager.hh"
#include "G4RunManagerFactory.hh"
#include "TestHelpers/IOManager.h"
#include "TestHelpers/Particle.h"
#include "TestHelpers/ParticleContainer.h"

class TestEnvironment : public ::testing::Environment
{
public:
  // The Location of the python script to plot the transport
  inline static const std::string PYTHON_SCRIPT =
      std::string(PYTHON_SCRIPT_DIR) + "plot_transport.py";
  // The location of the simplified geometry file
  inline static const std::string GEO_FILE_PATH =
      std::string(TEST_GEO_DIR) + "simplified_geo.gdml";
  // Boolean flag to enable verbose tracking
  inline static const bool ENABLE_VERBOSE_TRACKING = false;
  // The list of particles to be tested
  inline static const TestHelpers::ParticleContainer PARTICLE_CONTAINER = []()
  {
    TestHelpers::ParticleContainer container;
    container.addParticle(
        22, 0.2, 0.3 * CLHEP::GeV, 0.0, G4ThreeVector(0.0, 0.0, 0.0));
    container.addParticle(
        22, 4.5, 0.3 * CLHEP::GeV, 0.0, G4ThreeVector(0.0, 0.0, 0.0));
    container.addParticle(
        2212, 1, 0.2 * CLHEP::GeV, 0.0, G4ThreeVector(0.0, 0.0, 0.0));
    return container;
  }();
};

class TransportTests : public ::testing::TestWithParam<TestHelpers::Particle>
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
    runManager->SetUserInitialization(new DetectorConstruction);

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
  const TestHelpers::Particle& particle = GetParam();
  // Create the output directory of the test
  const std::string output_dir =
      TestHelpers::IOManager::create_test_output_dir();

  // Create action initialization to set up particle to simulate
  auto* actionInitialization = new ActionInitialization();

  // Tell the action initialization to set the particl
  actionInitialization->set_particle(&particle);
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
  ASSERT_TRUE(system(("python3 " + TestEnvironment::PYTHON_SCRIPT + " --input "
                      + transport_output_path + " --output " + plot_output_path)
                         .c_str())
              == 0);
}

INSTANTIATE_TEST_SUITE_P(
    ParticleTransportTests,
    TransportTests,
    ::testing::ValuesIn(TestEnvironment::PARTICLE_CONTAINER.get()),
    testing::PrintToStringParamName());