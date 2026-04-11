// Copyright (c) 2024 CERN for the benefit of the FastCaloSim project

#include "DetectorConstruction.hh"

#include <G4MagneticField.hh>
#include "G4UniformMagField.hh"

#include "FastSimModel.hh"
#include "G4Box.hh"
#include "G4FieldManager.hh"
#include "G4LogicalVolume.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4RegionStore.hh"
#include "G4TransportationManager.hh"
#include "G4Tubs.hh"
#include "TestHelpers/GeoPluginLoader.h"
#include "TestHelpers/MagFieldPlugin.h"


DetectorConstruction::DetectorConstruction()
    : G4VUserDetectorConstruction()
{
}

DetectorConstruction::~DetectorConstruction() {}

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  G4NistManager* nist = G4NistManager::Instance();

  // Define materials
  G4Material* worldMat = nist->FindOrBuildMaterial("G4_AIR");

  // World volume
  G4double world_size = 100 * CLHEP::m;
  G4Box* solidWorld =
      new G4Box("World", world_size / 2, world_size / 2, world_size / 2);
  G4LogicalVolume* world_log =
      new G4LogicalVolume(solidWorld, worldMat, "World");
  G4VPhysicalVolume* world_phys = new G4PVPlacement(
      0, G4ThreeVector(), world_log, "World", 0, false, 999, true);

  // Define the fast simulation trigger region
  // Particles that reach this region will be processed by the fast simulation
  // model
  auto detectorRegion = new G4Region(fFastSimTriggerRegionName);

  // Add the cylinders of the ID-Calo boundary to the fast simulation trigger
  // region
  TestHelpers::IDCaloBoundary boundary;
  std::map<std::string, TestHelpers::Cylinder> cylinders = boundary.getCylinders();
  for (const auto& [name, cyl] : cylinders) {
    add_cylinder(world_log, detectorRegion, cyl, name);
  }

  return world_phys;
}

void DetectorConstruction::ConstructSDandField()
{
  // Retrieve the static Geant4 field manager
  G4FieldManager* fieldMgr =
        G4TransportationManager::GetTransportationManager()->GetFieldManager();

  if (fUseAtlasField){
    G4cout<<"*************************************************"<<G4endl;
    G4cout<< "          Using ATLAS magnetic field"           << G4endl;
    G4cout<<"*************************************************"<<G4endl;

    // Load the ATLAS magnetic field
    GeoPluginLoader<MagFieldPlugin> loader;
    MagFieldPlugin* plugin = loader.load(ATLASMagneticFieldMapPluginLib);
    G4MagneticField* field = plugin->getField();

    if (field == nullptr) {
        G4cerr << "Failed to load magnetic field plugin" << G4endl;
        abort();
    }
    fField.Put(field);

    fieldMgr->SetDetectorField(dynamic_cast<G4MagneticField*>(fField.Get()));
    fieldMgr->CreateChordFinder(fField.Get());
  } else {

    G4cout<<"*************************************************"<<G4endl;
    G4cout<< "          Using uniform magnetic field"          << G4endl;
    G4cout<<"*************************************************"<<G4endl;

    auto uniformMagField     = new G4UniformMagField(G4ThreeVector(0.0, 0.0, fMagFieldStrength));
    fieldMgr->SetDetectorField(uniformMagField);
    fieldMgr->CreateChordFinder(uniformMagField);

  }


  // Create the fast simulation model
  auto triggerRegion =
      G4RegionStore::GetInstance()->GetRegion(fFastSimTriggerRegionName);
  auto fastSimModel = new FastSimModel("model", triggerRegion);
}

void DetectorConstruction::add_cylinder(G4LogicalVolume* world_log,
                                        G4Region* region,
                                        TestHelpers::Cylinder cyl,
                                        const std::string& name)
{
  // Compute the cylinder half-length
  const float halfLength = std::abs(cyl.zmax - cyl.zmin) / 2;

  // Create the solid volume
  G4Tubs* solid =
      new G4Tubs(name, cyl.rmin, cyl.rmax, halfLength, 0, 2 * CLHEP::pi);
  // Create the logical volume
  G4LogicalVolume* logic =
      new G4LogicalVolume(solid, world_log->GetMaterial(), name);

  // Place a cylinder in the positive z half space
  new G4PVPlacement(nullptr,
                    G4ThreeVector(0, 0, cyl.zmin + halfLength),
                    logic,
                    name + "::POS",
                    world_log,
                    false,
                    0,
                    true);
  // Place a cylinder in the negative z half space
  new G4PVPlacement(nullptr,
                    G4ThreeVector(0, 0, -(cyl.zmin + halfLength)),
                    logic,
                    name + "::NEG",
                    world_log,
                    false,
                    1,
                    true);

  // Add the cylinder to the (fast simulation) region
  region->AddRootLogicalVolume(logic);
}
