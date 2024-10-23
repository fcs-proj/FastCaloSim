// Copyright (c) 2024 CERN for the benefit of the FastCaloSim project

#ifndef DETECTORCONSTRUCTION_HH
#define DETECTORCONSTRUCTION_HH

#include "CLHEP/Units/SystemOfUnits.h"
#include "G4Cache.hh"
#include "G4VUserDetectorConstruction.hh"
#include "TestHelpers/IDCaloBoundary.h"

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4GlobalMagFieldMessenger;
class G4MagneticField;
class G4Region;


class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
  DetectorConstruction();
  virtual ~DetectorConstruction() override;

  virtual G4VPhysicalVolume* Construct() override;
  virtual void ConstructSDandField() override;

  // Interface to access the name of the fast simulation trigger region
  auto getFastSimTriggerRegionName() const -> std::string
  {
    return fFastSimTriggerRegionName;
  }
  // Set the use of the ATLAS magnetic field
  // If set to false, a uniform magnetic field is used
  void setUseAtlasField(bool useAtlasField)
  {
    fUseAtlasField = useAtlasField;
  }
  // Set the magnetic field strength of the uniform magnetic field
  // Only used if ATLAS magnetic field is not used
  void setMagFieldStrength(G4double magFieldStrength)
  {
    fMagFieldStrength = magFieldStrength;
  }

private:
  // Pointer to the thread-local fields
  G4Cache<G4MagneticField*> fField;
  // Name of the fast simulation trigger region
  const std::string fFastSimTriggerRegionName = "FastSimTriggerRegion";
  // Method to add a cylinder to the fast simulation trigger region
  void add_cylinder(G4LogicalVolume* world_log,
                    G4Region* region,
                    TestHelpers::Cylinder cyl,
                    const std::string& name);
  // Use ATLAS magnetic field
  bool fUseAtlasField = true;
  // Magnetic field strength of uniform magnetic field if ATLAS field is not used
  G4double fMagFieldStrength = 1 * CLHEP::tesla;
};

#endif
