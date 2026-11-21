// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#ifndef FASTSIMMODEL_HH
#define FASTSIMMODEL_HH

#include "FastCaloSim/Transport/G4CaloTransportTool.h"
#include "FastCaloSim/Extrapolation/FastCaloSimCaloExtrapolation.h"
#include "FastCaloSim/Core/TFCSSimulationState.h"
#include "FastCaloSim/Core/TFCSParametrizationBase.h"
#include <CLHEP/Random/RanluxEngine.h>


#include <utility>

#include <G4FieldTrack.hh>
#include <G4String.hh>
#include <G4ThreeVector.hh>
#include <G4Types.hh>
#include "G4VFastSimulationModel.hh"
#include "TestHelpers/TrackContainer.h"
#include "TestHelpers/ExtrapStateContainer.h"
#include "TestHelpers/SimStateContainer.h"

class G4FastTrack;
class G4ParticleDefinition;
class G4Region;
class G4VPhysicalVolume;
class G4FieldTrack;
class G4Navigator;

class FastSimModel : public G4VFastSimulationModel
{
 public:
  // The model transports particles through transportWorld, which is shared by
  // the models of all threads
  FastSimModel(const G4String& /*aModelName*/, G4Region* /*aEnvelope*/, G4VPhysicalVolume* transportWorld);
  ~FastSimModel() override;
  FastSimModel(const FastSimModel&) = delete;
  auto operator=(const FastSimModel&) -> FastSimModel& = delete;
  FastSimModel(FastSimModel&&) = delete;
  auto operator=(FastSimModel&&) -> FastSimModel& = delete;

  auto ModelTrigger(const G4FastTrack& /*unused*/) -> G4bool final;
  auto IsApplicable(const G4ParticleDefinition& /*unused*/) -> G4bool final;
  void DoIt(const G4FastTrack& /*aFastTrack*/ /*unused*/, G4FastStep& /*aFastStep*/ /*unused*/) final;
  void DoStep(G4FieldTrack&, G4Navigator*);

  // Get the current track container
  auto GetTransportTracks() const -> TestHelpers::TrackContainer { return fTransportTracks; }
  // Get extrapolation states
  auto GetExtrapolations() const -> TestHelpers::ExtrapStateContainer { return fExtrapolations; }
  // Get simulation states
  auto GetSimulationStates() const -> TestHelpers::SimStateContainer { return fSimulationStates; }

  // Reset the current state of the model
  void Reset() {
    fTransportTracks.clear();
    fExtrapolations.clear();
    fSimulationStates.clear();
    }

  // Set the boolean flag whether to do simulation
  void doSimulation(G4bool doSimulation) { fDoSimulation = doSimulation; }
  // Set the extrapolator tool
  void setExtrapolationTool(FastCaloSimCaloExtrapolation tool) { fExtrapolationTool = std::move(tool); }
  // Set debug flag
  void setDebug(bool debug) { m_debug = debug; }
  // Set FastCaloSim parametrization
  void setParametrization(TFCSParametrizationBase* param) { fParametrization = param; }

  private:

  // Geant4 transport tool of this model's thread
  G4CaloTransportTool fTransportTool;
  // Extrapolation tool
  FastCaloSimCaloExtrapolation fExtrapolationTool;
  // Current vector of tracks
  TestHelpers::TrackContainer fTransportTracks;
  // Current vector of extrapolations
  TestHelpers::ExtrapStateContainer fExtrapolations;
  // Current vector of simulation states
  TestHelpers::SimStateContainer fSimulationStates;
  // Core FastCaloSim API
  TFCSParametrizationBase* fParametrization{nullptr};

  // Boolean flag whether to do the simulation
  G4bool fDoSimulation{true};

  // Random engine
  CLHEP::RanluxEngine m_random_engine;
  // Boolean flag to enable debug output
  bool m_debug{false};

};
#endif /* FASTSIMMODEL_HH */
