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
#include "TestHelpers/ParticleSampler.h"

namespace
{

// Location of the simplified geometry
const std::string kTestGeoFile =
    std::string(TEST_GEO_DIR) + "simplified_geo.gdml";
// Location of the python script to plot the particle transport
const std::string kPythonScript =
    std::string(PYTHON_SCRIPT_DIR) + "plot_transport.py";
// Properties of particles to test in the transport
const std::vector<int> kParticleIds = {22, 11, -11};
const std::vector<float> kEtaValues = {-0.8, 0.8};
const float kKineticEnergy = 1 * CLHEP::GeV;
// Note: we generate at the origin but we only simulate the primary particles in
// the fast simulation DoIt method
const G4ThreeVector kPrimaryVertex(0, 0, 0);

}  // namespace

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
    fParser.Read(kTestGeoFile, false);

    // Construct the default run manager
    runManager =
        G4RunManagerFactory::CreateRunManager(G4RunManagerType::Serial);

    // Construct the detector
    runManager->SetUserInitialization(new DetectorConstruction);

    // Enable verbose tracking
    bool enableVerboseTracking = false;
    G4TrackingManager* trackingManager =
        G4EventManager::GetEventManager()->GetTrackingManager();
    trackingManager->SetVerboseLevel(static_cast<G4int>(enableVerboseTracking));

    // Configure fast simulation
    auto* fastSimulationPhysics = new G4FastSimulationPhysics();
    fastSimulationPhysics->BeVerbose();
    fastSimulationPhysics->ActivateFastSimulation("e-");
    fastSimulationPhysics->ActivateFastSimulation("e+");
    fastSimulationPhysics->ActivateFastSimulation("gamma");

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
  ASSERT_TRUE(system(("python3 " + kPythonScript + " --input "
                      + transport_output_path + " --output " + plot_output_path)
                         .c_str())
              == 0);
}

// Calls the ParticleTransportTest for each particle returned by the
// create_particles method
INSTANTIATE_TEST_SUITE_P(
    ParticleTransportTests,
    TransportTests,
    ::testing::ValuesIn(TestHelpers::ParticleSampler::generate_particles(
        kParticleIds, kEtaValues, kKineticEnergy, 0, kPrimaryVertex)),
    testing::PrintToStringParamName());
