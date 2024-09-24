/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

//
// BFieldVectorZR.h
//
// Magnetic field value (Bz,Br) stored in 2d map.
//
// Masahiro Morii, Harvard University
//
#ifndef BFIELDVECTORZR_H
#define BFIELDVECTORZR_H

class BFieldVectorZR
{
public:
  // constructor
  BFieldVectorZR() { ; }
  BFieldVectorZR(float Bz, float Br)
  {
    m_B[0] = Bz;
    m_B[1] = Br;
  }
  // setter
  void set(float Bz, float Br)
  {
    m_B[0] = Bz;
    m_B[1] = Br;
  }
  // accessors
  float z() const { return m_B[0]; }
  float r() const { return m_B[1]; }
  // array-like accessor
  float operator[](int i) const { return m_B[i]; }

private:
  float m_B[2];
};

#endif
