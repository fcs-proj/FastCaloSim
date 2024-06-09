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
  BFieldSolenoid()
      : m_orig(0)
      , m_tilt(0)
  {
    ;
  }
  // destructor
  ~BFieldSolenoid()
  {
    delete m_orig;
    if (m_orig != m_tilt)
      delete m_tilt;
  }
  // read/write map from/to file
  int readMap(std::istream& input);
  // int readMap( TFile* rootfile );
  // void writeMap( TFile* rootfile, bool tilted = false );
  //  move and tilt the map
  void moveMap(double dx, double dy, double dz, double ax, double ay);
  // compute magnetic field
  void getB(const double* xyz, double* B, double* deriv = 0) const;
  // accessor
  const BFieldMesh<double>* tiltedMap() const { return m_tilt; }

private:
  // data members
  BFieldMesh<double>* m_orig;  // original map as it was read from file
  BFieldMesh<double>* m_tilt;  // tilted and moved map
  // cache for speed
  mutable BFieldCache m_cache;
};

#endif
