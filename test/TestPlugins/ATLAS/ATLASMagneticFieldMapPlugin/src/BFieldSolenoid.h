/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

//
// BFieldSolenoid.h
//
// Magnetic field map for the ATLAS solenoid
//
// Masahiro Morii, Harvard University
//
#ifndef BFIELDSOLENOID_H
#define BFIELDSOLENOID_H

#include <iostream>
#include <vector>
// #include "TFile.h"
#include "BFieldZone.h"

class BFieldSolenoid
{
public:
  // constructor
  BFieldSolenoid() = default;
  // destructor
  ~BFieldSolenoid()
  {
    delete m_orig;
    if (m_orig != m_tilt) {
      delete m_tilt;
    }
  }
  // the maps are owned by this object
  BFieldSolenoid(const BFieldSolenoid&) = delete;
  auto operator=(const BFieldSolenoid&) -> BFieldSolenoid& = delete;
  BFieldSolenoid(BFieldSolenoid&&) = delete;
  auto operator=(BFieldSolenoid&&) -> BFieldSolenoid& = delete;
  // read/write map from/to file
  auto readMap(std::istream& input) -> int;
  // int readMap( TFile* rootfile );
  // void writeMap( TFile* rootfile, bool tilted = false );
  //  move and tilt the map
  void moveMap(double dx, double dy, double dz, double ax, double ay);
  // compute magnetic field
  void getB(const double* xyz, double* B, double* deriv = nullptr) const;
  // accessor
  auto tiltedMap() const -> const BFieldMesh<double>* { return m_tilt; }

private:
  // data members
  BFieldMesh<double>* m_orig {
      nullptr};  // original map as it was read from file
  BFieldMesh<double>* m_tilt {nullptr};  // tilted and moved map
  // cache for speed
  mutable BFieldCache m_cache;
};

#endif
