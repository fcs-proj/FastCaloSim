// Copyright (c) 2024 CERN for the benefit of the FastCaloSim project

#include "FastCaloSim/Transport/G4CaloTransportTool.h"

// Geant4 includes for for particle extrapolation
#include "G4FieldTrack.hh"
#include "G4FieldTrackUpdator.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4Navigator.hh"
#include "G4PVPlacement.hh"
#include "G4PathFinder.hh"
#include "G4TransportationManager.hh"

using namespace FastCaloSim::Transport;

G4CaloTransportTool::G4CaloTransportTool() = default;

G4CaloTransportTool::~G4CaloTransportTool()
{
  // Delete the world volume if we created it
  if (m_useSimplifiedGeo) {
    delete m_worldVolume;
  }

  // Delete the navigators and propagators for each thread
  for (auto& mapPair : m_propagatorHolder.getMap()) {
    delete mapPair.second->GetNavigatorForPropagating();
    delete mapPair.second;
  }
}

void G4CaloTransportTool::initializePropagator()
{
  G4cout << "Initializing G4PropagatorInField for thread "
         << G4Threading::G4GetThreadId() << G4endl;

  if (!m_worldVolume) {
    // If not set, get either the simplified or full world volume
    m_worldVolume = getWorldVolume();

    if (!m_worldVolume) {
      G4Exception("G4CaloTransportTool",
                  "FailedToGetWorldVolume",
                  FatalException,
                  "G4CaloTransportTool: Failed to get world volume.");
      abort();
    }
    G4cout << "Using world volume: " << m_worldVolume->GetName() << G4endl;
    G4cout << "Transport will be stopped at volume: " << m_transportLimitVolume
           << G4endl;
    G4cout << "Maximum allowed number of steps in particle transport: "
           << m_maxSteps << G4endl;
  }

  // Check if we already have propagator set up for the current thread
  auto* propagator = m_propagatorHolder.get();
  // If not, we create one
  if (!propagator) {
    propagator = makePropagator();
    m_propagatorHolder.set(propagator);
  } else {
    G4cerr << "G4CaloTransportTool::initializePropagator() Propagator already "
              "initialized!"
           << G4endl;
  }
}

auto G4CaloTransportTool::getWorldVolume() -> G4VPhysicalVolume*
{
  if (m_useSimplifiedGeo) {
    G4cout << "Creating simplified world volume for particle transport"
           << G4endl;
    // Get the logical world volume of the simplified geometry by name
    G4LogicalVolume* logVol = G4LogicalVolumeStore::GetInstance()->GetVolume(
        m_simplifiedWorldLogName);

    // Create the physical volume of the simplified world
    return new G4PVPlacement(
        nullptr,  // no rotation
        G4ThreeVector(0, 0, 0),  // world center at (0,0,0)
        logVol,  // logical volume
        "simplifiedWorldPhysVol",  // name of physical volume
        nullptr,  // mother volume
        false,  // not used
        999,  // copy number
        false);  // overlap check

  } else {
    G4cout << "Using full geometry for particle transport" << G4endl;
    return G4TransportationManager::GetTransportationManager()
        ->GetNavigatorForTracking()
        ->GetWorldVolume();
  }
}

auto G4CaloTransportTool::makePropagator() -> G4PropagatorInField*
{
  // Create a new navigator
  G4Navigator* navigator = new G4Navigator();
  // Set world volume in which the navigator will operate
  navigator->SetWorldVolume(m_worldVolume);
  // Get the global field manager
  G4FieldManager* fieldMgr =
      G4TransportationManager::GetTransportationManager()->GetFieldManager();
  // Create a new magnetic field propagator
  G4PropagatorInField* propagator =
      new G4PropagatorInField(navigator, fieldMgr);

  return propagator;
}

void G4CaloTransportTool::doStep(G4FieldTrack& fieldTrack)
{
  // Get the propagator and navigator for the current thread
  auto propagator = m_propagatorHolder.get();
  auto navigator = propagator->GetNavigatorForPropagating();

  G4double retSafety = -1.0;
  G4double currentMinimumStep = 10.0 * CLHEP::m;

  G4VPhysicalVolume* currentPhysVol =
      navigator->LocateGlobalPointAndSetup(fieldTrack.GetPosition(), nullptr);

  G4ThreeVector direction = fieldTrack.GetMomentumDirection();
  // Must be called before calling the computeStep method
  navigator->LocateGlobalPointAndSetup(fieldTrack.GetPosition(), &direction);

  if (fieldTrack.GetCharge() == 0) {
    /* Neutral particles: transport with navigator */

    // Compute the step length
    G4double stepLength =
        navigator->ComputeStep(fieldTrack.GetPosition(),
                               fieldTrack.GetMomentumDirection(),
                               currentMinimumStep,
                               retSafety);

    // Update the position of the track from the computed step length
    fieldTrack.SetPosition(fieldTrack.GetPosition()
                           + stepLength
                               * fieldTrack.GetMomentumDirection().unit());

  } else {
    /* Charged particles: transport with magnetic field propagator */
    propagator->ComputeStep(
        fieldTrack, currentMinimumStep, retSafety, currentPhysVol);
  }
}

std::vector<G4FieldTrack> G4CaloTransportTool::transport(
    const G4Track& G4InputTrack)
{
  // Get the navigator for the current thread
  auto* navigator = m_propagatorHolder.get()->GetNavigatorForPropagating();

  // Create a vector to store the output steps
  std::vector<G4FieldTrack> outputStepVector;
  // Initialize the tmpFieldTrack with the input track
  G4FieldTrack tmpFieldTrack('0');
  G4FieldTrackUpdator::Update(&tmpFieldTrack, &G4InputTrack);
  // Fill with the initial particle position
  outputStepVector.push_back(tmpFieldTrack);

  // Iterate until we reach the maximum number of steps or the requested volume
  for (unsigned int iStep = 0; iStep < m_maxSteps; iStep++) {
    // Perform a single Geant4 step
    doStep(tmpFieldTrack);
    // Fill the output vector with the updated track
    outputStepVector.push_back(tmpFieldTrack);
    // Get the name of the volume in which the particle is located
    auto volume = navigator->LocateGlobalPointAndSetup(
        tmpFieldTrack.GetPosition(), nullptr);
    if (volume != nullptr) {
      std::string volName = volume->GetName();
      // We stop the track navigation once we have reached the provided volume
      if (volName.find(m_transportLimitVolume) != std::string::npos) {
        break;
      }
    } else {
      G4ExceptionDescription description;
      description
          << "Transport failure at step " << iStep << "/" << m_maxSteps
          << G4endl << " - Position: " << tmpFieldTrack.GetPosition() << G4endl
          << " - Momentum: " << tmpFieldTrack.GetMomentum() << G4endl
          << "Possible cause: The transport is likely outside the world volume."
          << G4endl
          << "Check if an envelope volume is defined and properly set up."
          << G4endl << "This issue should not occur during normal operation.";

      G4Exception("G4CaloTransportTool::transport",
                  "LocateGlobalPointAndSetup failed: Particle may be "
                  "transported outside the world volume.",
                  JustWarning,
                  description);
      break;
    }
  }

  return outputStepVector;
}
