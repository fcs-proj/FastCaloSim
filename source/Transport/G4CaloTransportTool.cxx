// Copyright (c) 2024 CERN for the benefit of the FastCaloSim project

#include "FastCaloSim/Transport/G4CaloTransportTool.h"

// Geant4 includes for for particle extrapolation
#include <iomanip>
#include <iostream>
#include <sstream>

#include "FastCaloSim/Core/FCSDebugPrint.h"
#include "G4AtlasRK4.hh"
#include "G4ChordFinder.hh"
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
#include "G4VIntegrationDriver.hh"

namespace
{
// RAII guard: forces G4AtlasRK4 to recompute its cached momentum-to-force
// conversion factor (m_cof/m_imom, derived from the momentum magnitude
// m_mom) on every Stepper() call for the lifetime of this guard, then
// restores whatever tolerance was configured before. See the long comment
// in transport() for why this exists. No-op if the stepper isn't a
// G4AtlasRK4 (e.g. a differently-configured job) or is null.
struct ScopedAtlasRK4CacheDisable
{
  G4AtlasRK4* stepper;
  G4double savedFraction;

  explicit ScopedAtlasRK4CacheDisable(G4AtlasRK4* s)
      : stepper(s)
      , savedFraction(s ? s->GetChangeFraction() : 0.0)
  {
    if (stepper) {
      stepper->SetChangeFraction(0.0);
    }
  }

  ~ScopedAtlasRK4CacheDisable()
  {
    if (stepper) {
      stepper->SetChangeFraction(savedFraction);
    }
  }

  ScopedAtlasRK4CacheDisable(const ScopedAtlasRK4CacheDisable&) = delete;
  ScopedAtlasRK4CacheDisable& operator=(const ScopedAtlasRK4CacheDisable&) =
      delete;
};
}  // namespace

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
  // threads (as the previous per-worker initializePropagator() did) is a data
  // race that can produce duplicate world volumes and subtly different
  // navigation between threads. Building it once on the master thread also
  // keeps the Geant4 split-class per-thread data for the world correctly sized
  // across all workers. std::call_once guarantees once-only execution but not
  // which thread runs it, so restrict the actual creation to the master thread;
  // worker threads only report whether the shared world is ready.
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
  // initializeGeometry() before any worker transports a particle. transport()
  // now builds its own navigator/propagator per call (see the rationale there),
  // so there is nothing to cache here; this is purely a readiness check.
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
  // Rationale: previously one propagator/navigator was cached per thread and
  // reset at the top of every call. Making that reset bit-for-bit equivalent to
  // a freshly constructed propagator proved impractical -- Geant4 spreads
  // per-track state across the propagator (looping/zero-step counters,
  // intersection-locator caches, first/last-step-in-volume flags), its
  // navigator (located-point and safety caches, zero-step/push machinery) and
  // the field manager's chord finder, and any field the reset misses leaks
  // between calls. Because the mapping of events to worker threads differs
  // between single- and multi-threaded running, such leakage makes a shower
  // depend on which events the thread processed before it, which breaks single-
  // vs multi-threaded reproducibility. Constructing the state fresh per call
  // removes that dependency by construction: each transport() sees only the
  // input track. The cost (one navigator + propagator per transported particle)
  // is negligible next to the field-integration steps below, and it mirrors how
  // every Geant4 worker owns a private navigator over the shared, read-only
  // geometry. The objects are stack-scoped and destruct in reverse declaration
  // order (propagator before navigator), so there is nothing to clean up.
  G4Navigator navigator;
  navigator.SetWorldVolume(m_worldVolume);
  G4FieldManager* fieldMgr =
      G4TransportationManager::GetTransportationManager()->GetFieldManager();
  G4PropagatorInField propagator(&navigator, fieldMgr);

  // The chord finder lives on the (thread-local) field manager, not on the
  // freshly built propagator, so it is the one piece of state not already reset
  // by construction above. Clear its cached last step-size estimate so field
  // propagation does not depend on the previous transport on this thread.
  G4FieldManagerStore::GetInstance()->ClearAllChordFindersState();

  // The stepper itself (owned by the chord finder, which lives on the
  // thread-persistent field manager -- NOT rebuilt fresh above, and NOT
  // touched by ClearAllChordFindersState()) can carry its own leftover
  // per-track state across calls. Concretely: G4AtlasRK4 (ATLAS's default
  // field stepper, "Sim.G4Stepper"="AtlasRK4") caches the momentum
  // magnitude of whichever track last used it (m_mom, and the derived
  // m_cof/m_imom force-scaling factors) and only recomputes them if the new
  // track's momentum differs by more than a relative tolerance
  // (GetChangeFraction(), default 1e-7) -- a valid optimization *within* one
  // track's own step sequence (momentum magnitude barely changes step to
  // step in a magnetic field) but not across the boundary into an unrelated
  // track. Since transport() is called once per particle and the stepper is
  // shared across the whole thread's lifetime, if this new particle's
  // momentum happens to coincidentally fall within that tolerance of
  // whatever the thread's previous particle left behind, the stale
  // force-scaling factor gets reused for this particle's first RK4
  // sub-steps -- a tiny perturbation that an adaptive stepper responds to
  // with a different step count, without changing where the track
  // ultimately ends up (matches the observed symptom exactly: identical
  // IDCaloBoundary crossing point, different caloSteps.size() between ST
  // and MT for a small fraction of particles).
  G4AtlasRK4* atlasStepper = nullptr;
  if (G4ChordFinder* chordFinder = fieldMgr->GetChordFinder()) {
    if (G4VIntegrationDriver* driver = chordFinder->GetIntegrationDriver()) {
      atlasStepper = dynamic_cast<G4AtlasRK4*>(driver->GetStepper());
    }
  }

  // FCS_DEBUG_EXTRAPOL_EVENT + FCS_DEBUG_POS_FILTER_FILE-gated confirmation:
  // log the current particle's momentum against whatever is cached in the
  // stepper *before* the fix below touches anything, to directly show
  // whether the two coincidentally nearly match (the trigger condition for
  // the bug described above).
  if (atlasStepper && std::getenv("FCS_DEBUG_EXTRAPOL_EVENT")
      && fcsDebugShouldPrintPos(G4InputTrack.GetPosition().x(),
                                G4InputTrack.GetPosition().y(),
                                G4InputTrack.GetPosition().z()))
  {
    std::ostringstream oss;
    oss << std::setprecision(17) << "[FCS_DEBUG_STEPPER_CACHE] joinPos=("
        << G4InputTrack.GetPosition().x() << ","
        << G4InputTrack.GetPosition().y() << ","
        << G4InputTrack.GetPosition().z()
        << ") curMom=" << G4InputTrack.GetMomentum().mag()
        << " cachedMom=" << atlasStepper->GetCachedMomentum()
        << " fraction=" << atlasStepper->GetChangeFraction() << "\n";
    std::lock_guard<std::mutex> lock(fcsDebugPrintMutex());
    std::cout << oss.str() << std::flush;
  }

  // The fix: force every Stepper() call during this transport() to recompute
  // the momentum-derived scaling factors fresh from this particle's own
  // momentum, regardless of what the stepper's cache holds on entry. This
  // makes this transport() call's result depend only on its own input, never
  // on which particle happened to run on this thread before it. Restored to
  // whatever it was on scope exit, so this only affects FastCaloSim's own
  // re-transport, not the shared field manager's behavior for any other
  // (e.g. real G4 tracking) use on this thread.
  //
  // FCS_DISABLE_STEPPER_FIX=1 skips applying it (confirmation print above
  // still fires), so the natural bug can be caught directly instead of
  // simultaneously being masked by this same run's fix.
  ScopedAtlasRK4CacheDisable stepperCacheGuard(
      std::getenv("FCS_DISABLE_STEPPER_FIX") ? nullptr : atlasStepper);

  // Initialize the tmpFieldTrack with the input track
  G4FieldTrack tmpFieldTrack('0');
  G4FieldTrackUpdator::Update(&tmpFieldTrack, &G4InputTrack);

  // FCS_DEBUG_EXTRAPOL_EVENT + FCS_DEBUG_POS_FILTER_FILE-gated: dump every
  // step's position for specific, already-identified particles (matched by
  // their pre-transport position) to localize exactly where an ST-vs-MT
  // caloSteps.size() mismatch (same input, same final IDCaloBoundary
  // crossing, different step count) first diverges.
  const bool fcsDebugThisTransport = std::getenv("FCS_DEBUG_EXTRAPOL_EVENT")
      && fcsDebugShouldPrintPos(G4InputTrack.GetPosition().x(),
                                G4InputTrack.GetPosition().y(),
                                G4InputTrack.GetPosition().z());
  const auto fcsDebugPrintStep = [&](int stepIdx)
  {
    if (!fcsDebugThisTransport)
      return;
    std::ostringstream oss;
    oss << std::setprecision(17) << "[FCS_DEBUG_TRANSPORT] joinPos=("
        << G4InputTrack.GetPosition().x() << ","
        << G4InputTrack.GetPosition().y() << ","
        << G4InputTrack.GetPosition().z() << ") step=" << stepIdx << " pos=("
        << tmpFieldTrack.GetPosition().x() << ","
        << tmpFieldTrack.GetPosition().y() << ","
        << tmpFieldTrack.GetPosition().z() << ")\n";
    std::lock_guard<std::mutex> lock(fcsDebugPrintMutex());
    std::cout << oss.str() << std::flush;
  };

  // Establish the starting location with an absolute (non-relative) search so
  // the first step does not depend on any prior navigator history.
  navigator.LocateGlobalPointAndSetup(
      tmpFieldTrack.GetPosition(), nullptr, /*relativeSearch=*/false);

  // Fill with the initial particle position
  outputStepVector.push_back(tmpFieldTrack);
  fcsDebugPrintStep(0);

  // Iterate until we reach the maximum number of steps or the requested volume
  for (unsigned int iStep = 0; iStep < m_maxSteps; iStep++) {
    // Perform a single Geant4 step
    doStep(propagator, tmpFieldTrack);
    // Fill the output vector with the updated track
    outputStepVector.push_back(tmpFieldTrack);
    fcsDebugPrintStep(static_cast<int>(iStep) + 1);
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
