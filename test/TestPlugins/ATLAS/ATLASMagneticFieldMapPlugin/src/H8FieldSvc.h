/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

//
// H8FieldSvc.h
//
// Magnetic field map of the H8 beam test
//
// Masahiro Morii, Harvard University
//
#ifndef H8FIELDSVC_H
#define H8FIELDSVC_H

// FrameWork includes
// #include "AthenaBaseComps/AthService.h"

// MagField includes
#include "BFieldH8Grid.h"
#include "IMagFieldSvc.h"

// STL includes
#include <string>
#include <vector>

// forward declarations
class CondAttrListCollection;

namespace MagField
{
// class H8FieldSvc{ : public extends<AthService, IMagFieldSvc> {
class H8FieldSvc : public IMagFieldSvc
{
public:
  //** Constructor with parameters */
  explicit H8FieldSvc(const std::string& name /*, ISvcLocator* pSvcLocator*/);

  /** Destructor */
  ~H8FieldSvc() override;
  H8FieldSvc(const H8FieldSvc&) = delete;
  auto operator=(const H8FieldSvc&) -> H8FieldSvc& = delete;
  H8FieldSvc(H8FieldSvc&&) = delete;
  auto operator=(H8FieldSvc&&) -> H8FieldSvc& = delete;

  /** Athena algorithm's interface methods */
  static auto initialize() -> bool;
  auto start() -> bool;
  static auto finalize() -> bool;

  /** get B field value at given position */
  /** xyz[3] is in mm, bxyz[3] is in kT */
  /** if deriv[9] is given, field derivatives are returned in kT/mm */
  void getField(const double* xyz,
                double* B,
                double* deriv = nullptr) const override;
  /** getFieldZR simply calls getField **/
  void getFieldZR(const double* xyz,
                  double* B,
                  double* deriv = nullptr) const override;

private:
  // initialize map
  auto initializeMap() -> bool;
  // read the field map
  auto readMap(const std::string& mapFile) -> bool;

  /** Data members **/

  // field map name
  std::string m_H8MapFilename;

  // field data
  std::vector<BFieldH8Grid> m_grid;

  // first magnet offset
  double m_dx1 {0};
  double m_dy1 {0};
  double m_dz1 {0};
  // second magnet offset
  double m_dx2 {43830};
  double m_dy2 {0};
  double m_dz2 {0};
};
}  // namespace MagField

#endif
