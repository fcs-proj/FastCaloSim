/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

#include <iostream>

#include "AtlasFieldSvc.h"
#include "G4MagneticField.hh"
#include "IMagFieldSvc.h"
#include "MagFieldPlugin.h"

namespace
{
class AtlasField : public G4MagneticField
{
public:
  // Construct the field object from the IMagFieldSvc
  explicit AtlasField(MagField::IMagFieldSvc* mfield);

  // Implementation of G4 method to retrieve field value
  void GetFieldValue(const double* point, double* field) const override
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
            << m_magFieldSvc_AtlasField << " to " << mfield << '\n';
}

class ATLASMagneticFieldMapPlugin : public MagFieldPlugin

{
public:
  // Constructor
  ATLASMagneticFieldMapPlugin();

  // Destructor
  ~ATLASMagneticFieldMapPlugin() override;

  ATLASMagneticFieldMapPlugin(const ATLASMagneticFieldMapPlugin&) = delete;
  auto operator=(const ATLASMagneticFieldMapPlugin&)
      -> ATLASMagneticFieldMapPlugin& = delete;
  ATLASMagneticFieldMapPlugin(ATLASMagneticFieldMapPlugin&&) = delete;
  auto operator=(ATLASMagneticFieldMapPlugin&&)
      -> ATLASMagneticFieldMapPlugin& = delete;

  // Overriding virtual function
  auto getField() -> G4MagneticField* override;

private:
  MagField::IMagFieldSvc* m_magFieldSvc {nullptr};
};

ATLASMagneticFieldMapPlugin::ATLASMagneticFieldMapPlugin()
{
  std::cout << "HELLO from Atlas Magnetic Field Plugin" << '\n';
}

ATLASMagneticFieldMapPlugin::~ATLASMagneticFieldMapPlugin()
{
  std::cout << "GOODBYE from Atlas Magnetic Field Plugin" << '\n';
}

auto ATLASMagneticFieldMapPlugin::getField() -> G4MagneticField*
{
  auto* atlasFieldSvs = new MagField::AtlasFieldSvc(
      TEST_BASE_DIR
          + std::string("data/fieldmap/bmagatlas_09_fullAsym20400.data"),
      /*isAscii=*/true);
  atlasFieldSvs->handle();
  m_magFieldSvc = atlasFieldSvs;

  std::cout << "AtlasFieldSvc::makeField with m_magFieldSvc " << m_magFieldSvc
            << '\n';
  return new AtlasField(&*m_magFieldSvc);
}
}  // namespace

extern "C" auto createATLASMagneticFieldMapPlugin() -> MagFieldPlugin*
{
  return new ATLASMagneticFieldMapPlugin;
}
