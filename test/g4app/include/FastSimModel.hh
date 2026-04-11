// Copyright (c) 2024 CERN for the benefit of the FastCaloSim project

#ifndef FASTSIMMODEL_HH
#define FASTSIMMODEL_HH

#include "FastCaloSim/Transport/G4CaloTransportTool.h"
#include "FastCaloSim/Extrapolation/FastCaloSimCaloExtrapolation.h"
#include "FastCaloSim/Core/TFCSSimulationState.h"
#include "FastCaloSim/Core/TFCSParametrizationBase.h"
#include <CLHEP/Random/RanluxEngine.h>


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
class G4FieldTrack;
class G4Navigator;

using namespace FastCaloSim::Core;
using FastCaloSim::Extrapolation::FastCaloSimCaloExtrapolation;
using FastCaloSim::Transport::G4CaloTransportTool;


class FastSimModel : public G4VFastSimulationModel
{
 public:
  FastSimModel(G4String, G4Region*);
  ~FastSimModel();

  virtual G4bool ModelTrigger(const G4FastTrack&) final;
  virtual G4bool IsApplicable(const G4ParticleDefinition&) final;
  virtual void DoIt(const G4FastTrack&, G4FastStep&) final;
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
  void setExtrapolationTool(FastCaloSimCaloExtrapolation tool) { fExtrapolationTool = tool; }
  // Set debug flag
  void setDebug(bool debug) { m_debug = debug; }
  // Set FastCaloSim parametrization
  void setParametrization(TFCSParametrizationBase* param) { fParametrization = param; }

  private:

  // Geant4 transport tool
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
  TFCSParametrizationBase* fParametrization;

  // Boolean flag whether to do the simulation
  G4bool fDoSimulation;

  // Random engine
  CLHEP::RanluxEngine m_random_engine;
  // Boolean flag to enable debug output
  bool m_debug;

};
#endif /* FASTSIMMODEL_HH */
