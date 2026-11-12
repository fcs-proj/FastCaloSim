// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#include "FastCaloSim/Transport/G4CaloTransportTool.h"

// Geant4 includes for particle extrapolation
#include "G4ChargeState.hh"
#include "G4DynamicParticle.hh"
#include "G4EquationOfMotion.hh"
#include "G4FieldManagerStore.hh"
#include "G4FieldTrack.hh"
#include "G4FieldTrackUpdator.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4Navigator.hh"
#include "G4PVPlacement.hh"
#include "G4PropagatorInField.hh"
#include "G4Track.hh"
#include "G4TransportationManager.hh"

G4CaloTransportTool::G4CaloTransportTool(G4VPhysicalVolume* world)
{
  if (!world) {
    G4Exception("G4CaloTransportTool::G4CaloTransportTool",
                "MissingWorldVolume",
                FatalException,
                "G4CaloTransportTool: no world volume to transport in.");
  }
  // The navigator and propagator are built once and reset for every particle
  // in transport(): building them per particle took about an eighth of the
  // transport time.
  m_navigator = std::make_unique<G4Navigator>();
  m_navigator->SetWorldVolume(world);
  m_propagator = std::make_unique<G4PropagatorInField>(
      m_navigator.get(),
      G4TransportationManager::GetTransportationManager()->GetFieldManager());
}

G4CaloTransportTool::~G4CaloTransportTool() = default;

auto G4CaloTransportTool::createSimplifiedWorld(
    const std::string& worldLogicalVolumeName) -> G4VPhysicalVolume*
{
  G4LogicalVolume* logVol =
      G4LogicalVolumeStore::GetInstance()->GetVolume(worldLogicalVolumeName);
  if (!logVol) {
    G4ExceptionDescription description;
    description << "G4CaloTransportTool: simplified world logical volume '"
                << worldLogicalVolumeName
                << "' was not found in the G4LogicalVolumeStore. Ensure the "
                   "simplified transport geometry is loaded first.";
    G4Exception("G4CaloTransportTool::createSimplifiedWorld",
                "MissingSimplifiedWorldLog",
                FatalException,
                description);
  }

  // Create the physical volume of the simplified world. Creating it on the
  // master thread keeps the Geant4 split-class per-thread data of the world
  // correctly sized across all workers.
  G4cout << "Creating simplified world volume for particle transport" << G4endl;
  return new G4PVPlacement(nullptr,  // no rotation
                           G4ThreeVector(0, 0, 0),  // world center at (0,0,0)
                           logVol,  // logical volume
                           "simplifiedWorldPhysVol",  // name of physical volume
                           nullptr,  // mother volume
                           false,  // not used
                           999,  // copy number
                           false);  // overlap check
}

void G4CaloTransportTool::doStep(G4FieldTrack& fieldTrack,
                                 G4VPhysicalVolume* currentPhysVol)
{
  G4double retSafety = -1.0;
  G4double currentMinimumStep = 10.0 * CLHEP::m;

  if (fieldTrack.GetCharge() == 0) {
    /* Neutral particles: transport with navigator */

    // Compute the step length
    G4double stepLength =
        m_navigator->ComputeStep(fieldTrack.GetPosition(),
                                 fieldTrack.GetMomentumDirection(),
                                 currentMinimumStep,
                                 retSafety);

    // Update the position of the track from the computed step length
    fieldTrack.SetPosition(fieldTrack.GetPosition()
                           + stepLength
                               * fieldTrack.GetMomentumDirection().unit());

  } else {
    /* Charged particles: transport with magnetic field propagator */
    m_propagator->ComputeStep(
        fieldTrack, currentMinimumStep, retSafety, currentPhysVol);
  }
}

auto G4CaloTransportTool::transport(const G4Track& G4InputTrack)
    -> std::vector<G4FieldTrack>
{
  // Create a vector to store the output steps
  std::vector<G4FieldTrack> outputStepVector;

  // Reset the navigator and propagator to the state of freshly built ones.
  // The propagator uses the thread's current global field manager, both as
  // the detector one and as the current one (which provides the equation of
  // motion below), as a freshly built propagator does: a framework may replace
  // the field manager after this tool was built (Athena does when it sets up
  // its magnetic field). ClearPropagatorState() resets the propagator's
  // per-track state (looping and zero-step counters, first-step flags),
  // ResetStackAndState() the navigator's history and state, and the chord
  // finder is cleared further down. Each transport() call therefore sees only
  // its own input track, whatever was transported before, which keeps single-
  // and multi-threaded results identical.
  m_propagator->SetDetectorFieldManager(
      G4TransportationManager::GetTransportationManager()->GetFieldManager());
  m_propagator->FindAndSetFieldManager(nullptr);
  m_propagator->ClearPropagatorState();
  m_navigator->ResetStackAndState();

  // Configure the (thread-persistent, shared) equation of motion with THIS
  // particle's charge/mass/moments before stepping, exactly as real Geant4
  // tracking does on every step (see G4Transportation::AlongStepDoIt, which
  // builds a fresh G4ChargeState from the current track and calls
  // SetChargeMomentumMass() every time). The equation of motion is shared
  // across the thread's entire lifetime and is not reset by
  // ClearPropagatorState() above, so skipping this would leave it holding
  // whatever charge/mass the previous track on this thread left behind (real
  // G4 tracking of some other particle, or FastCaloSim's own previous
  // transport() call). A wrong charge flips the Lorentz-force curvature
  // outright -- a completely different trajectory from the very first step,
  // not a subtle numerical difference -- so this must run before every
  // transport() call regardless of what ran before it.
  const G4DynamicParticle* inputDynamicParticle =
      G4InputTrack.GetDynamicParticle();
  const G4ParticleDefinition* inputParticleDef =
      inputDynamicParticle->GetDefinition();
  G4ChargeState chargeState(inputDynamicParticle->GetCharge(),
                            inputDynamicParticle->GetMagneticMoment(),
                            inputParticleDef->GetPDGSpin());
  m_propagator->GetCurrentEquationOfMotion()->SetChargeMomentumMass(
      chargeState,
      inputDynamicParticle->GetTotalMomentum(),
      inputDynamicParticle->GetMass());

  // The chord finder lives on the (thread-local) field manager, not on the
  // propagator. Clear its cached last step-size estimate so field propagation
  // does not depend on the previous transport on this thread.
  G4FieldManagerStore::GetInstance()->ClearAllChordFindersState();

  // Initialize the tmpFieldTrack with the input track
  G4FieldTrack tmpFieldTrack('0');
  G4FieldTrackUpdator::Update(&tmpFieldTrack, &G4InputTrack);

  // Locate the starting volume with an absolute (non-relative) search
  G4VPhysicalVolume* volume = m_navigator->LocateGlobalPointAndSetup(
      tmpFieldTrack.GetPosition(), nullptr, /*relativeSearch=*/false);

  // Fill with the initial particle position
  outputStepVector.push_back(tmpFieldTrack);

  // Iterate until we reach the maximum number of steps or the requested volume
  for (unsigned int iStep = 0; iStep < m_maxSteps; iStep++) {
    // Perform a single Geant4 step from the volume located last
    doStep(tmpFieldTrack, volume);
    // Fill the output vector with the updated track
    outputStepVector.push_back(tmpFieldTrack);
    // Locate the volume the particle is now in, which the next step starts
    // from. The direction resolves the volume if the step ended on an edge.
    const G4ThreeVector direction = tmpFieldTrack.GetMomentumDirection();
    volume = m_navigator->LocateGlobalPointAndSetup(tmpFieldTrack.GetPosition(),
                                                    &direction);
    if (volume != nullptr) {
      // We stop the track navigation once we have reached the provided volume
      if (volume->GetName().find(m_transportLimitVolume) != std::string::npos) {
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
