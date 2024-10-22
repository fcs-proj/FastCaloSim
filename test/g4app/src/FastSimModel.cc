// Copyright (c) 2024 CERN for the benefit of the FastCaloSim project

#include "FastSimModel.hh"
#include <G4FastStep.hh>
#include <G4FastTrack.hh>
#include <G4Track.hh>
#include <G4VFastSimulationModel.hh>
#include <G4Region.hh>
#include "G4FieldTrackUpdator.hh"
#include "G4FieldTrack.hh"
#include "G4RunManager.hh"
#include "FastCaloSim/Core/TFCSTruthState.h"
#include "FastCaloSim/Core/TFCSExtrapolationState.h"


class G4ParticleDefinition;


FastSimModel::FastSimModel(G4String aModelName, G4Region* aEnvelope)
  : G4VFastSimulationModel(aModelName, aEnvelope)
  , fDoSimulation(true)
  , m_debug(false)
  , fParametrization(nullptr)
{
  fTransportTool.initializePropagator();
  m_random_engine.setSeed(42);
}

FastSimModel::~FastSimModel() {}


G4bool FastSimModel::IsApplicable(const G4ParticleDefinition& aParticleType)
{
  return true;
}


G4bool FastSimModel::ModelTrigger(const G4FastTrack& aFastTrack)
{
  return true;
}



void FastSimModel::DoIt(const G4FastTrack& aFastTrack, G4FastStep& aFastStep)
{

  // Get Geant4 primary track
  const G4Track * track = aFastTrack.GetPrimaryTrack();

  // Set the FastCaloSim truth state
  TFCSTruthState truth;
  truth.set_pdgid( track -> GetDefinition() -> GetPDGEncoding());

  // Set the kinematics of the FastCaloSim truth state
  truth.SetPtEtaPhiM(track -> GetMomentum().perp(),
                    track -> GetMomentum().eta(),
                    track -> GetMomentum().phi(),
                      track -> GetDefinition() -> GetPDGMass());

  // Set the vertex of the FastCaloSim truth state
  truth.set_vertex( track -> GetPosition().x(),
                    track -> GetPosition().y(),
                    track -> GetPosition().z());

  if(m_debug){
    std::cout << "[FastSimModel::DoIt] Received particle with pid="
              << truth.pdgid() << std::endl;

    std::cout << "[FastSimModel::DoIt] Particle has position x=" << truth.X()
              << " y=" << truth.Y()
              << " z=" << truth.Z()
              << " eta=" << truth.Eta()
              << " phi=" << truth.Phi()
              << " r=" << truth.Perp()
              << " ekin=" << truth.Ekin() << std::endl;

    std::cout << "[FastSimModel::DoIt] Particle has momentum px=" << track->GetMomentum().x()
              << " py=" << track->GetMomentum().y()
              << " pz=" << track->GetMomentum().z()
              << " p=" << track->GetMomentum().mag()
              << " pt=" << track->GetMomentum().perp()
              << " eta=" << track->GetMomentum().eta() << std::endl;
  }

  // Perform the FastCaloSim transportation
  std::vector<G4FieldTrack> step_vector = fTransportTool.transport(*aFastTrack.GetPrimaryTrack());
  // Add the track to the vector of tracks
  TestHelpers::Track trk( step_vector );
  fTransportTracks.add(trk);

  // Do the extrapolation
  TFCSExtrapolationState extrap;
  fExtrapolationTool.extrapolate(extrap, &truth, step_vector);
  fExtrapolations.add(extrap);

  // Do the simulation if requested
  if (fDoSimulation) {

    if(!fParametrization){
      std::cerr << "[FastSimModel::DoIt] No parametrization set!" << std::endl;
      return;
    }
    // set the random engine
    TFCSSimulationState simul;
    simul.setRandomEngine(&m_random_engine);

    // Simulate the energy response of the particle
    fParametrization->simulate(simul, &truth, &extrap);

    // Add the simulation state to the vector of simulation states
    fSimulationStates.add(simul);
  }

  // Kill particle
  aFastStep.KillPrimaryTrack();
  aFastStep.SetPrimaryTrackPathLength(0.0);

}
