#include <iostream>

#include "AtlasFieldSvc.h"
#include "G4MagneticField.hh"
#include "IMagFieldSvc.h"
#include "MagFieldPlugin.h"

class AtlasField : public G4MagneticField
{
public:
  // Construct the field object from the IMagFieldSvc
  AtlasField(MagField::IMagFieldSvc* m);

  // Implementation of G4 method to retrieve field value
  void GetFieldValue(const double* point, double* field) const
  {
    m_magFieldSvc_AtlasField->getField(point, field);
  }

private:
  // Pointer to the magnetic field service.
  // We use a raw pointer here to avoid ServiceHandle overhead.
  MagField::IMagFieldSvc* m_magFieldSvc_AtlasField;
};

AtlasField::AtlasField(MagField::IMagFieldSvc* mfield)
    : m_magFieldSvc_AtlasField(mfield)
{
  std::cout << "New instance of AtlasField, setting m_magFieldSvc_AtlasField: "
            << m_magFieldSvc_AtlasField << " to " << mfield << std::endl;
}

class ATLASMagneticFieldMapPlugin : public MagFieldPlugin

{
public:
  // Constructor
  ATLASMagneticFieldMapPlugin();

  // Destructor
  ~ATLASMagneticFieldMapPlugin();

protected:
  // Overriding virtual function
  G4MagneticField* getField();

private:
  MagField::IMagFieldSvc* m_magFieldSvc;
};

ATLASMagneticFieldMapPlugin::ATLASMagneticFieldMapPlugin()
{
  std::cout << "HELLO from Atlas Magnetic Field Plugin" << std::endl;
}

ATLASMagneticFieldMapPlugin::~ATLASMagneticFieldMapPlugin()
{
  std::cout << "GOODBYE from Atlas Magnetic Field Plugin" << std::endl;
}

G4MagneticField* ATLASMagneticFieldMapPlugin::getField()
{
  MagField::AtlasFieldSvc* atlasFieldSvs = new MagField::AtlasFieldSvc(
      (TEST_MAG_FIELD_DIR + std::string("bmagatlas_09_fullAsym20400.data"))
          .c_str(),
      true);
  atlasFieldSvs->handle();
  m_magFieldSvc = atlasFieldSvs;

  std::cout << "AtlasFieldSvc::makeField with m_magFieldSvc " << m_magFieldSvc
            << std::endl;
  return new AtlasField(&*m_magFieldSvc);
}

extern "C" ATLASMagneticFieldMapPlugin* createATLASMagneticFieldMapPlugin()
{
  return new ATLASMagneticFieldMapPlugin;
}
