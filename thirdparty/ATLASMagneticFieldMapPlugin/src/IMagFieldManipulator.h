/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

///////////////////////////////////////////////////////////////////
// IMagFieldManipulator.h, (c) ATLAS Detector software
///////////////////////////////////////////////////////////////////
#ifndef MAGFIELDINTERFACES_IMAGFIELDMANIPULATOR_H
#define MAGFIELDINTERFACES_IMAGFIELDMANIPULATOR_H

#include <cmath>
#include <iostream>

// Framework includes
// #include "GaudiKernel/IAlgTool.h"

// Amg classes
// #include "GeoPrimitives/GeoPrimitives.h"

/** Declaration of the interface ID ( interface id, major version, minor
 * version) */
// static const InterfaceID IID_IMagFieldManipulator("IMagFieldManipulator", 1,
// 0);

namespace MagField
{

/** @ class IMagFieldManipulator

 @ author Valerio.Ippolito -at- cern.ch
 */
// class IMagFieldManipulator: virtual public IAlgTool {

class IMagFieldManipulator
{
  ///////////////////////////////////////////////////////////////////
  // Public methods:
  ///////////////////////////////////////////////////////////////////
public:
  /** constructor */
  IMagFieldManipulator() { ; }

  /** Retrieve interface ID */
  //        static const InterfaceID& interfaceID() {
  //            return IID_IMagFieldManipulator;
  //        }

  /** change the point where the field should be evaluated */
  /** in AtlasFieldSvc, this is called before B field is evaluated and */
  /** fed into modifyField **/
  /** xyz_new[3] and xyz_old[3] are in mm */
  virtual void modifyPosition(const double* xyz_old, double* xyz_new) = 0;

  /** correct B field value at a position xyz */
  /** bxyz[3] is in kT */
  /** if deriv[9] is given, field derivatives are returned in kT/mm */
  virtual void modifyField(double* bxyz, double* deriv = 0) = 0;
};
}  // namespace MagField

#endif  //> !MAGFIELDINTERFACES_IMAGFIELDMANIPULATOR_H
