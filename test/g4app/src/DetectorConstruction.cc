#include "G4SystemOfUnits.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4RegionStore.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "FastSimModel.hh"
#include "DetectorConstruction.hh"
#include "G4Tubs.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"
#include "G4GDMLParser.hh"


G4ThreadLocal
G4GlobalMagFieldMessenger* DetectorConstruction::fMagFieldMessenger = nullptr;

DetectorConstruction::DetectorConstruction()
    : G4VUserDetectorConstruction()
{
}

DetectorConstruction::~DetectorConstruction()
{
}

G4VPhysicalVolume *DetectorConstruction::Construct()
{

  G4NistManager* nist = G4NistManager::Instance();
  
  // Define materials
  G4Material* worldMat = nist->FindOrBuildMaterial("G4_AIR");

  // World volume
  G4double world_size = 100 * m; 
  G4Box* solidWorld = new G4Box("World", world_size / 2, world_size / 2, world_size / 2);
  G4LogicalVolume* logicWorld = new G4LogicalVolume(solidWorld, worldMat, "World");
  G4VPhysicalVolume* physWorld = new G4PVPlacement(0, G4ThreeVector(), logicWorld, "World", 0, false, 999, true);

  // Barrel cylinder definition
  const float cylinderThickness = 10*CLHEP::mm;
  const float barrelRMin = 1148*CLHEP::mm;
  const float barrelRMax = barrelRMin + cylinderThickness;
  const float barrelZMin = 0.0*CLHEP::mm;
  const float barrelZMax = 3549.5*CLHEP::mm;

  G4Tubs* solid = new G4Tubs("IDCaloBoundary::Barrel", barrelRMin, barrelRMax, std::abs(barrelZMax - barrelZMin), 0, 2 * CLHEP::pi);
  G4LogicalVolume* logic = new G4LogicalVolume(solid, worldMat, "IDCaloBoundary::Barrel");
  new G4PVPlacement(nullptr, G4ThreeVector(0, 0, 0), logic, "IDCaloBoundary::Barrel", logicWorld, false, 0, true);

  // Define the fast simulation trigger region
  // Particles that reach this region will be processed by the fast simulation model
  auto detectorRegion = new G4Region(fFastSimTriggerRegionName);
  // Add the barrel cylinder to the region
  detectorRegion->AddRootLogicalVolume(logic);

  return physWorld;
}

void DetectorConstruction::ConstructSDandField()
{ 
  // Set a magnetic field of 1 Tesla in x direction
  fMagFieldMessenger = new G4GlobalMagFieldMessenger(G4ThreeVector(0, 1*CLHEP::tesla, 0));

  // Register the field messenger for deleting
  G4AutoDelete::Register(fMagFieldMessenger);

  // Create the fast simulation model
  auto triggerRegion = G4RegionStore::GetInstance()->GetRegion(fFastSimTriggerRegionName);
  auto fastSimModel = new FastSimModel("model", triggerRegion);

}