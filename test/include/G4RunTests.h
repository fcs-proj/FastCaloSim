#pragma once

#include <gtest/gtest.h>

#include "ActionInitialization.hh"
#include "DetectorConstruction.hh"
#include "FTFP_BERT.hh"
#include "G4FastSimulationPhysics.hh"
#include "G4GDMLParser.hh"
#include "G4RunManagerFactory.hh"
#include "TestConfig/G4RunTestsConfig.h"
#include "TestHelpers/Event.h"

class G4RunTests : public virtual ::testing::TestWithParam<TestHelpers::Event>
{
protected:
  // Make the run manager a member variable so that it can be accessed in the
  // tests
  static G4RunManager* run_manager;
  static ActionInitialization* actionInitialization;

  // Per-test-suite set-up.
  // Note: the geant4 kernel does not support kernel re-instanciation
  // So this needs to be done once per test suite
  static void SetUpTestSuite()
  {
    // Load the simplified geometry from a GDML file into the volume store
    G4GDMLParser fParser;
    fParser.Read(G4RunTestsConfig::GEO_FILE_PATH, false);

    // Construct the default run manager
    run_manager =
        G4RunManagerFactory::CreateRunManager(G4RunManagerType::Serial);

    // Construct the detector
    auto* detectorConstruction = new DetectorConstruction();
    detectorConstruction->setUseAtlasField(G4RunTestsConfig::USE_ATLAS_FIELD);
    run_manager->SetUserInitialization(detectorConstruction);

    // Set verbosity of the tracking manager
    G4TrackingManager* trackingManager =
        G4EventManager::GetEventManager()->GetTrackingManager();
    trackingManager->SetVerboseLevel(
        static_cast<G4int>(G4RunTestsConfig::ENABLE_VERBOSE_TRACKING));

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
    run_manager->SetUserInitialization(physicsList);

    // Create action initialization to set up event to simulate
    actionInitialization = new ActionInitialization();
  }

  // Per-test-suite tear-down
  static void TearDownTestSuite()
  {
    // Clean up the run manager
    delete run_manager;
  }
};

// Initialize the static members
G4RunManager* G4RunTests::run_manager = nullptr;
ActionInitialization* G4RunTests::actionInitialization = nullptr;