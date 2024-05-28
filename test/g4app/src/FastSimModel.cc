//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
#include "FastSimModel.hh"
#include <G4FastStep.hh>                 
#include <G4FastTrack.hh>                
#include <G4Track.hh>                    
#include <G4VFastSimulationModel.hh>     
#include <G4Region.hh>                
#include "FastCaloSim/Transport/G4CaloTransportTool.h" 
#include "G4FieldTrackUpdator.hh"         
#include "G4FieldTrack.hh"                


class G4ParticleDefinition;


FastSimModel::FastSimModel(G4String aModelName, G4Region* aEnvelope)
  : G4VFastSimulationModel(aModelName, aEnvelope)
{}

FastSimModel::FastSimModel(G4String aModelName)
  : G4VFastSimulationModel(aModelName)
{}


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
  // To make our life easier, we only consider primary particles in tests
  if (aFastTrack.GetPrimaryTrack()->GetParentID() != 0) {
    aFastStep.KillPrimaryTrack();
    aFastStep.SetPrimaryTrackPathLength(0.0);
    return;
  }
  float eta = aFastTrack.GetPrimaryTrack()->GetMomentum().eta();
  float phi = aFastTrack.GetPrimaryTrack()->GetMomentum().phi();
  float r = aFastTrack.GetPrimaryTrack()->GetPosition().perp();
  float ekin = aFastTrack.GetPrimaryTrack()->GetKineticEnergy();

  std::cout<<"[FastSimModel::DoIt] Received particle with eta="<<eta<< " phi="<<phi<<" r="<<r<<" ekin="<<ekin<<std::endl;

  G4CaloTransportTool transportTool;
  transportTool.initializePropagator();

  std::vector<G4FieldTrack> step_vector = transportTool.transport(*aFastTrack.GetPrimaryTrack());

  // Save the step vector as csv in case that saveTransport was set
  if (!fTransportOutputPath.empty()) {
    std::ofstream file;
    file.open(fTransportOutputPath);
    file << "x,y,z" << std::endl;
    for (auto& step : step_vector) {
      file << step.GetPosition().x() << "," << step.GetPosition().y() << "," << step.GetPosition().z() << std::endl;
    }
    file.close();
  }

  // Kill particle
  aFastStep.KillPrimaryTrack();
  aFastStep.SetPrimaryTrackPathLength(0.0);

}