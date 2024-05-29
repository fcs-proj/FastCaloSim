#ifndef FASTSIMMODEL_HH
#define FASTSIMMODEL_HH

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

  void saveTransport(std::string outputPath) { fTransportOutputPath = outputPath; }
  std::string fTransportOutputPath = "";

};
#endif /* FASTSIMMODEL_HH */
