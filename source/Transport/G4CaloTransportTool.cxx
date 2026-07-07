// Copyright (c) 2024 CERN for the benefit of the FastCaloSim project

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
#include "G4Threading.hh"
#include "G4Track.hh"
#include "G4TransportationManager.hh"

G4CaloTransportTool::G4CaloTransportTool() = default;

G4CaloTransportTool::~G4CaloTransportTool()
{
  // Delete the world volume if we created it. The per-call navigator and
  // propagator built inside transport() are stack-scoped and clean up
  // themselves, so there is no per-thread transport state to release here.
  if (m_useSimplifiedGeo) {
    delete m_worldVolume;
  }
}

bool G4CaloTransportTool::initializeGeometry()
{
  // The world volume is shared by all threads' navigators and must be created
  // exactly once, on the master thread. getWorldVolume() registers a
  // G4PVPlacement into the global Geant4 geometry stores in the
  // simplified-geometry case; doing so concurrently from multiple worker
  // threads is a data race that can produce duplicate world volumes and
  // subtly different navigation between threads. Building it once on the
  // master thread also keeps the Geant4 split-class per-thread data for the
  // world correctly sized across all workers. std::call_once guarantees
  // once-only execution but not which thread runs it, so restrict the actual
  // creation to the master thread; worker threads only report whether the
  // shared world is ready.
  if (!G4Threading::IsMasterThread()) {
    return m_worldVolume != nullptr;
  }

  std::call_once(m_worldVolumeOnceFlag,
                 [this]()
                 {
                   // Get either the simplified or full world volume
                   m_worldVolume = getWorldVolume();

                   if (!m_worldVolume) {
                     G4Exception("G4CaloTransportTool::initializeGeometry",
                                 "FailedToGetWorldVolume",
                                 JustWarning,
                                 "G4CaloTransportTool: Failed to get world "
                                 "volume.");
                     return;
                   }
                   G4cout << "Using world volume: " << m_worldVolume->GetName()
                          << G4endl;
                   G4cout << "Transport will be stopped at volume: "
                          << m_transportLimitVolume << G4endl;
                   G4cout << "Maximum allowed number of steps in particle "
                             "transport: "
                          << m_maxSteps << G4endl;
                 });

  return m_worldVolume != nullptr;
}

bool G4CaloTransportTool::initializePropagator()
{
  // The shared world volume must have been created on the master thread by
  // initializeGeometry() before any worker transports a particle. Each call
  // to transport() builds its own navigator/propagator (see the comment
  // there), so this is purely a readiness check.
  if (!m_worldVolume) {
    G4Exception("G4CaloTransportTool::initializePropagator",
                "WorldVolumeNotInitialized",
                JustWarning,
                "G4CaloTransportTool: world volume is not initialized. "
                "initializeGeometry() must be called on the master thread "
                "before initializePropagator().");
    return false;
  }

  return true;
}

auto G4CaloTransportTool::getWorldVolume() -> G4VPhysicalVolume*
{
  if (m_useSimplifiedGeo) {
    G4cout << "Creating simplified world volume for particle transport"
           << G4endl;
    // Get the logical world volume of the simplified geometry by name
    G4LogicalVolume* logVol = G4LogicalVolumeStore::GetInstance()->GetVolume(
        m_simplifiedWorldLogName);

    if (!logVol) {
      G4ExceptionDescription description;
      description << "G4CaloTransportTool: simplified world logical volume '"
                  << m_simplifiedWorldLogName
                  << "' was not found in the G4LogicalVolumeStore. Ensure the "
                     "simplified transport geometry is loaded before "
                     "initializeGeometry() runs.";
      G4Exception("G4CaloTransportTool::getWorldVolume",
                  "MissingSimplifiedWorldLog",
                  JustWarning,
                  description);
      return nullptr;
    }

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

void G4CaloTransportTool::doStep(G4PropagatorInField& propagator,
                                 G4FieldTrack& fieldTrack)
{
  // Use the navigator owned by the propagator built for this transport() call
  auto* navigator = propagator.GetNavigatorForPropagating();

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
    propagator.ComputeStep(
        fieldTrack, currentMinimumStep, retSafety, currentPhysVol);
  }
}

std::vector<G4FieldTrack> G4CaloTransportTool::transport(
    const G4Track& G4InputTrack)
{
  // Create a vector to store the output steps
  std::vector<G4FieldTrack> outputStepVector;

  // The shared world volume must have been created on the master thread by
  // initializeGeometry() before any particle is transported. Fail loudly
  // instead of navigating a null world.
  if (!m_worldVolume) {
    G4Exception("G4CaloTransportTool::transport",
                "WorldVolumeNotInitialized",
                JustWarning,
                "G4CaloTransportTool: world volume is not initialized. "
                "initializeGeometry() must be called on the master thread "
                "before transport().");
    return outputStepVector;
  }

  // Build a fresh navigator and propagator for this single transport() call.
  //
  // Geant4 spreads per-track state across the propagator (looping/zero-step
  // counters, intersection-locator caches, first/last-step-in-volume flags),
  // its navigator (located-point and safety caches, zero-step/push
  // machinery), and the field manager's chord finder. A propagator/navigator
  // reused across calls would carry any of this state forward unless every
  // field were explicitly reset, and a missed field leaks state between
  // calls. Since the mapping of events to worker threads differs between
  // single- and multi-threaded running, such leakage would make a shower's
  // result depend on which events the thread happened to process before it --
  // breaking single- vs multi-threaded reproducibility. Constructing the
  // state fresh per call removes that dependency by construction: each
  // transport() call sees only its own input track. The cost (one navigator
  // + propagator per transported particle) is negligible next to the
  // field-integration steps below, and it mirrors how every Geant4 worker
  // owns a private navigator over the shared, read-only geometry. The
  // objects are stack-scoped and destruct in reverse declaration order
  // (propagator before navigator), so there is nothing to clean up.
  G4Navigator navigator;
  navigator.SetWorldVolume(m_worldVolume);
  G4FieldManager* fieldMgr =
      G4TransportationManager::GetTransportationManager()->GetFieldManager();
  G4PropagatorInField propagator(&navigator, fieldMgr);

  // Configure the (thread-persistent, shared) equation of motion with THIS
  // particle's charge/mass/moments before stepping, exactly as real Geant4
  // tracking does on every step (see G4Transportation::AlongStepDoIt, which
  // builds a fresh G4ChargeState from the current track and calls
  // SetChargeMomentumMass() every time). The equation of motion is shared
  // across the thread's entire lifetime and is not rebuilt by the fresh
  // navigator/propagator constructed above, so skipping this would leave it
  // holding whatever charge/mass the previous track on this thread left
  // behind (real G4 tracking of some other particle, or FastCaloSim's own
  // previous transport() call). A wrong charge flips the Lorentz-force
  // curvature outright -- a completely different trajectory from the very
  // first step, not a subtle numerical difference -- so this must run before
  // every transport() call regardless of what ran before it.
  const G4DynamicParticle* inputDynamicParticle =
      G4InputTrack.GetDynamicParticle();
  const G4ParticleDefinition* inputParticleDef =
      inputDynamicParticle->GetDefinition();
  G4ChargeState chargeState(inputDynamicParticle->GetCharge(),
                            inputDynamicParticle->GetMagneticMoment(),
                            inputParticleDef->GetPDGSpin());
  propagator.GetCurrentEquationOfMotion()->SetChargeMomentumMass(
      chargeState,
      inputDynamicParticle->GetTotalMomentum(),
      inputDynamicParticle->GetMass());

  // The chord finder lives on the (thread-local) field manager, not on the
  // freshly built propagator, so it is the one piece of state not already reset
  // by construction above. Clear its cached last step-size estimate so field
  // propagation does not depend on the previous transport on this thread.
  G4FieldManagerStore::GetInstance()->ClearAllChordFindersState();

  // Initialize the tmpFieldTrack with the input track
  G4FieldTrack tmpFieldTrack('0');
  G4FieldTrackUpdator::Update(&tmpFieldTrack, &G4InputTrack);

  // Establish the starting location with an absolute (non-relative) search so
  // the first step does not depend on any prior navigator history.
  navigator.LocateGlobalPointAndSetup(
      tmpFieldTrack.GetPosition(), nullptr, /*relativeSearch=*/false);

  // Fill with the initial particle position
  outputStepVector.push_back(tmpFieldTrack);

  // Iterate until we reach the maximum number of steps or the requested volume
  for (unsigned int iStep = 0; iStep < m_maxSteps; iStep++) {
    // Perform a single Geant4 step
    doStep(propagator, tmpFieldTrack);
    // Fill the output vector with the updated track
    outputStepVector.push_back(tmpFieldTrack);
    // Get the name of the volume in which the particle is located
    auto volume = navigator.LocateGlobalPointAndSetup(
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
