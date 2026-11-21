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

#include <array>
#include <cstddef>

class BFieldVectorZR
{
public:
  // constructor
  BFieldVectorZR() = default;
  BFieldVectorZR(float Bz, float Br)
      : m_B {Bz, Br}
  {
  }
  // setter
  void set(float Bz, float Br)
  {
    m_B[0] = Bz;
    m_B[1] = Br;
  }
  // accessors
  auto z() const -> float { return m_B[0]; }
  auto r() const -> float { return m_B[1]; }
  // array-like accessor
  auto operator[](std::size_t i) const -> float { return m_B[i]; }

private:
  std::array<float, 2> m_B {};
};

#endif
