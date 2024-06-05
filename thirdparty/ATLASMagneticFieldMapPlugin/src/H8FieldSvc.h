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
  H8FieldSvc(const std::string& name /*, ISvcLocator* pSvcLocator*/);

  /** Destructor */
  virtual ~H8FieldSvc();

  /** Athena algorithm's interface methods */
  bool initialize();
  bool start();
  bool finalize();

  /** get B field value at given position */
  /** xyz[3] is in mm, bxyz[3] is in kT */
  /** if deriv[9] is given, field derivatives are returned in kT/mm */
  virtual void getField(const double* xyz,
                        double* bxyz,
                        double* deriv = nullptr) const;
  /** getFieldZR simply calls getField **/
  virtual void getFieldZR(const double* xyz,
                          double* bxyz,
                          double* deriv = nullptr) const;

private:
  // initialize map
  bool initializeMap();
  // read the field map
  bool readMap(const std::string mapFile);

  /** Data members **/

  // field map name
  std::string m_H8MapFilename;

  // field data
  std::vector<BFieldH8Grid> m_grid;

  // first magnet offset
  double m_dx1;
  double m_dy1;
  double m_dz1;
  // second magnet offset
  double m_dx2;
  double m_dy2;
  double m_dz2;
};
}  // namespace MagField

#endif
