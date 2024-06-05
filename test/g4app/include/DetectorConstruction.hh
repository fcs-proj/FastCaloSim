#ifndef DETECTORCONSTRUCTION_HH
#define DETECTORCONSTRUCTION_HH

#include "CLHEP/Units/SystemOfUnits.h"
#include "G4Cache.hh"
#include "G4VUserDetectorConstruction.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4GlobalMagFieldMessenger;
class G4MagneticField;
class G4Region;

struct Cylinder
{
  double rmin;
  double rmax;
  double zmin;
  double zmax;
};

struct IDCaloBoundaryDimension
{
  /**
   * @struct IdCaloBoundaryDimension
   * @brief Dimensions of a virtual boundary between Inner Detector and
   * calorimeter system of the ATLAS detector
   * @details The boundary is defined by a series of cylinders with 10 mm of
   * thickness that cover the ATLAS ID-Calo boundary. The values describe the
   * positive z half sapce of the boundary. The full boundary is constructed by
   * mirroring the positive z half space in the negative z half space.
   */
  // Thickness of the boundary
  const float cylinderThickness = 10 * CLHEP::mm;

  // Barrel cylinder definition
  Cylinder barrel = {1148 * CLHEP::mm,
                     1148 * CLHEP::mm + cylinderThickness,
                     0.0 * CLHEP::mm,
                     3550.0 * CLHEP::mm};
  // Endcap plates barrel
  Cylinder barrelEndcap = {120 * CLHEP::mm,
                           barrel.rmax,
                           barrel.zmax,
                           barrel.zmax + cylinderThickness};
  // Inner beampipe cylinder definition
  Cylinder innerBeamPipe = {barrelEndcap.rmin,
                            barrelEndcap.rmin + cylinderThickness,
                            barrelEndcap.zmax,
                            4587.0 * CLHEP::mm};
  // Endcap plates beam pipe
  Cylinder innerBeamPipeEndcap = {41 * CLHEP::mm,
                                  innerBeamPipe.rmin,
                                  innerBeamPipe.zmax,
                                  innerBeamPipe.zmax + cylinderThickness};
  // Outer beampipe cylinder definition
  Cylinder outerBeamPipe = {innerBeamPipeEndcap.rmin,
                            innerBeamPipeEndcap.rmin + cylinderThickness,
                            innerBeamPipeEndcap.zmax,
                            6783 * CLHEP::mm};
};

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

private:
  // Pointer to the thread-local fields
  G4Cache<G4MagneticField*> fField;
  // Name of the fast simulation trigger region
  const std::string fFastSimTriggerRegionName = "FastSimTriggerRegion";
  // Method to add a cylinder to the fast simulation trigger region
  void add_cylinder(G4LogicalVolume* world_log,
                    G4Region* region,
                    Cylinder cyl,
                    const std::string& name);
};

#endif
