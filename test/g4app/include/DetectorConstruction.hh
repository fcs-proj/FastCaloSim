#ifndef DETECTORCONSTRUCTION_HH
#define DETECTORCONSTRUCTION_HH

#include "G4VUserDetectorConstruction.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4GlobalMagFieldMessenger;

using CylinderDefinitionMap = std::map<std::string, std::tuple<float, float, float, float>>;

class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
    DetectorConstruction();
    virtual ~DetectorConstruction() override;

    virtual G4VPhysicalVolume *Construct() override;
    virtual void ConstructSDandField() override;

    auto getFastSimTriggerRegionName() const -> std::string { return fFastSimTriggerRegionName; }
    
private:
    // Magnetic field messenger
    static G4ThreadLocal G4GlobalMagFieldMessenger*  fMagFieldMessenger;
    const std::string fFastSimTriggerRegionName = "FastSimTriggerRegion";
};

#endif
