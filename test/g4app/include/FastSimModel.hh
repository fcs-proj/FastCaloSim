#ifndef FASTSIMMODEL_HH
#define FASTSIMMODEL_HH

#include "FastCaloSim/Transport/G4CaloTransportTool.h"


#include <G4FieldTrack.hh>
#include <G4String.hh>                // for G4String
#include <G4ThreeVector.hh>           // for G4ThreeVector
#include <G4Types.hh>                 // for G4bool, G4double
#include "G4VFastSimulationModel.hh"  // for G4VFastSimulationModel

class G4FastTrack;
class G4ParticleDefinition;
class G4Region;
class G4FieldTrack;
class G4Navigator;

class FastSimModel : public G4VFastSimulationModel
{
 public:
  FastSimModel(G4String, G4Region*);
  FastSimModel(G4String);
  ~FastSimModel();

  virtual G4bool ModelTrigger(const G4FastTrack&) final;
  virtual G4bool IsApplicable(const G4ParticleDefinition&) final;
  virtual void DoIt(const G4FastTrack&, G4FastStep&) final;
  void DoStep(G4FieldTrack&, G4Navigator*);

  void saveTransport(std::string outputPath) { 
    fTransportOutputPath = outputPath;
    fNEventMap[fTransportOutputPath] = 0;
    fNFastTrackMap[fTransportOutputPath] = 0;
    }
    
  std::string fTransportOutputPath = "";

  private:
    // Current event number for a specific transport path
    // As we are re-using the same class instance for multiple tests, we map the event number to the transport path
    std::map<std::string, int> fNEventMap;
    // Current number of tracks processed by this model for the current event and transport path
    std::map<std::string, int> fNFastTrackMap;
    // Writes the transport vector to a csv file
    void writeTransportData(std::vector<G4FieldTrack> step_vector);
    // Geant4 transport tool
    G4CaloTransportTool fTransportTool;

};
#endif /* FASTSIMMODEL_HH */
