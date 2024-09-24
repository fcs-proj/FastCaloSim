/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

///////////////////////////////////////////////////////////////////
// IMagFieldSvc.h, (c) ATLAS Detector software
///////////////////////////////////////////////////////////////////
#ifndef MAGFIELDINTERFACES_IMAGFIELDSVC_H
#define MAGFIELDINTERFACES_IMAGFIELDSVC_H

#include <cmath>
#include <iostream>

// Framework includes
// #include "GaudiKernel/IInterface.h"

// Amg classes
// #include "GeoPrimitives/GeoPrimitives.h"

namespace MagField
{

/** @ class IMagFieldSvc

 @ author Elmar.Ritsch -at- cern.ch
 */
// class IMagFieldSvc: virtual public IInterface {
class IMagFieldSvc
{
  ///////////////////////////////////////////////////////////////////
  // Public methods:
  ///////////////////////////////////////////////////////////////////
public:
  /** Creates the InterfaceID and interfaceID() method */
  // DeclareInterfaceID(IMagFieldSvc, 1, 0);

  /** constructor */
  IMagFieldSvc()
      : m_solenoidCurrent(0.0)
      , m_toroidsCurrent(0.0)
  {
    ;
  }

  /** get B field value at given position */
  /** xyz[3] is in mm, bxyz[3] is in kT */
  /** if deriv[9] is given, field derivatives are returned in kT/mm */
  virtual void getField(const double* xyz,
                        double* bxyz,
                        double* deriv = nullptr) const = 0;

  /** a getField() wrapper for Amg classes */
  // void getField(const Amg::Vector3D *xyz, Amg::Vector3D *bxyz) const {
  //   getField( xyz->data(), bxyz->data(), nullptr );
  // }
  // void getField(const Amg::Vector3D *xyz, Amg::Vector3D *bxyz,
  // Amg::RotationMatrix3D *deriv) const {
  //   getField( xyz->data(), bxyz->data(), deriv->data() );
  // }

  /** get B field value on the z-r plane at given position */
  /** works only inside the solenoid; otherwise calls getField() above */
  /** xyz[3] is in mm, bxyz[3] is in kT */
  /** if deriv[9] is given, field derivatives are returned in kT/mm */
  virtual void getFieldZR(const double* xyz,
                          double* bxyz,
                          double* deriv = nullptr) const = 0;

  /** status of the magnets */
  bool solenoidOn() const { return getSolenoidCurrent() > 0.0; }
  bool toroidsOn() const { return getToroidsCurrent() > 0.0; }
  float getSolenoidCurrent() const { return m_solenoidCurrent; }
  float getToroidsCurrent() const { return m_toroidsCurrent; }

  // protected:
  void setSolenoidCurrent(float current) { m_solenoidCurrent = current; }
  void setToroidsCurrent(float current) { m_toroidsCurrent = current; }

private:
  float m_solenoidCurrent;  // solenoid current in ampere
  float m_toroidsCurrent;  // toroid current in ampere
};
}  // namespace MagField

#endif  //> !MAGFIELDINTERFACES_IMAGFIELDSVC_H
