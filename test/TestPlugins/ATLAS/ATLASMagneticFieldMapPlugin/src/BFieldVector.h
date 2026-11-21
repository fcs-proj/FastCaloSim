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

#include <array>
#include <cstddef>

template<class T>
class BFieldVector
{
public:
  // constructor
  BFieldVector() = default;
  BFieldVector(T Bz, T Br, T Bphi)
      : m_B {Bz, Br, Bphi}
  {
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
  auto operator[](std::size_t i) const -> T { return m_B[i]; }

private:
  std::array<T, 3> m_B {};
};

#endif
