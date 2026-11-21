/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

//
// BFieldVector.h
//
// Magnetic field value stored in the map.
// It may be either short (toroid) or double (solenoid).
//
// Masahiro Morii, Harvard University
//
#ifndef BFIELDVECTOR_H
#define BFIELDVECTOR_H

template<class T>
class BFieldVector
{
public:
  // constructor
  BFieldVector() { ; }
  BFieldVector(T Bz, T Br, T Bphi)
  {
    m_B[0] = Bz;
    m_B[1] = Br;
    m_B[2] = Bphi;
  }
  // setter
  void set(T Bz, T Br, T Bphi)
  {
    m_B[0] = Bz;
    m_B[1] = Br;
    m_B[2] = Bphi;
  }
  // accessors
  auto z() const -> T { return m_B[0]; }
  auto r() const -> T { return m_B[1]; }
  auto phi() const -> T { return m_B[2]; }
  // array-like accessor
  auto operator[](int i) const -> T { return m_B[i]; }

private:
  T m_B[3];
};

#endif
