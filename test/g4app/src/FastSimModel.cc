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
#include "G4RunManager.hh"

class G4ParticleDefinition;


FastSimModel::FastSimModel(G4String aModelName, G4Region* aEnvelope)
  : G4VFastSimulationModel(aModelName, aEnvelope)
{
}

FastSimModel::FastSimModel(G4String aModelName)
  : G4VFastSimulationModel(aModelName)
{
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
  // Current event number
  G4int n_event = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
  
  // Counts the number of tracks processed by this model for the current event and test iteration
  if(n_event != fNEventMap[fTransportOutputPath]){
    fNEventMap[fTransportOutputPath] = n_event;
    fNFastTrackMap[fTransportOutputPath] = 0;
  }

  int pid = aFastTrack.GetPrimaryTrack()->GetParticleDefinition()->GetPDGEncoding();
  float eta = aFastTrack.GetPrimaryTrack()->GetPosition().eta();
  float phi = aFastTrack.GetPrimaryTrack()->GetPosition().phi();
  float r = aFastTrack.GetPrimaryTrack()->GetPosition().perp();
  float ekin = aFastTrack.GetPrimaryTrack()->GetKineticEnergy();

  std::cout<<"[FastSimModel::DoIt] Processing event "<<n_event<<" track "<<fNFastTrackMap[fTransportOutputPath]<<std::endl;
  std::cout<<"[FastSimModel::DoIt] Received particle with pid=" << pid << " eta="<<eta<< " phi="<<phi<<" r="<<r<<" ekin="<<ekin<<std::endl;

  G4CaloTransportTool transportTool;
  transportTool.initializePropagator();

  std::vector<G4FieldTrack> step_vector = transportTool.transport(*aFastTrack.GetPrimaryTrack());

  // write the transport data to a file if requested
  if (!fTransportOutputPath.empty()) { 
    writeTransportData(step_vector);
  }

  // Kill particle
  aFastStep.KillPrimaryTrack();
  aFastStep.SetPrimaryTrackPathLength(0.0);

  // Increase processed track number
  fNFastTrackMap[fTransportOutputPath]++;

}

void FastSimModel::writeTransportData(std::vector<G4FieldTrack> step_vector){

  // Check if this is the first call to DoIt
  bool initCall = (fNEventMap[fTransportOutputPath] == 0 && fNFastTrackMap[fTransportOutputPath] == 0);

  // Save the step vector as csv in case that saveTransport was set
  std::ofstream transport_file;
  if (initCall) {
    // If its the first call, create a new file
    transport_file.open(fTransportOutputPath);
  } else {
    // If its not the first call, append to the existing file
    transport_file.open(fTransportOutputPath, std::ios::app);
  } 

  if (!transport_file.is_open()) {
    G4cerr << "Error: Unable to open transport output file: " << fTransportOutputPath;
  }
  
  if(fNEventMap[fTransportOutputPath] == 0 && fNFastTrackMap[fTransportOutputPath] == 0){
    transport_file << "event,track,x,y,z,r,eta,phi" << std::endl;
  }

  for (auto& step : step_vector) {
    transport_file  << fNEventMap[fTransportOutputPath] << "," 
                    << fNFastTrackMap[fTransportOutputPath] << ","
                    << step.GetPosition().x() << "," 
                    << step.GetPosition().y() << "," 
                    << step.GetPosition().z() << ","
                    << step.GetPosition().perp() << ","
                    << step.GetPosition().eta() << ","
                    << step.GetPosition().phi()
                    << std::endl;
  }
  transport_file.close();
  
}