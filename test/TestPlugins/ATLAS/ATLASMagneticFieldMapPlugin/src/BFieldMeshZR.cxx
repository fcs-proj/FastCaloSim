/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

//
// BFieldMeshZR.cxx
//
// A 2-dim z-r mesh inside the solenoid field map
//
// Masahiro Morii, Harvard University
//
#include "BFieldMeshZR.h"

//
// Construct the look-up table to accelerate bin-finding.
//
void BFieldMeshZR::buildLUT()
{
  for (int j = 0; j < 2; j++) {  // z, r
    // determine the unit size, q, to be used in the LUTs
    double const width = m_mesh[j].back() - m_mesh[j].front();
    double q(width);
    for (unsigned i = 0; i < m_mesh[j].size() - 1; i++) {
      q = std::min(q, m_mesh[j][i + 1] - m_mesh[j][i]);
    }
    // find the number of units in the LUT
    int n = int(width / q) + 1;
    q = width / (n + 0.5);
    m_invUnit[j] = 1.0 / q;  // new unit size
    n++;
    int m = 0;  // mesh number
    m_LUT[j].reserve(static_cast<std::size_t>(n));
    for (int i = 0; i < n; i++) {  // LUT index
      if (i * q + m_mesh[j].front()
          > m_mesh[j][static_cast<std::size_t>(m + 1)])
      {
        m++;
      }
      m_LUT[j].push_back(m);
    }
  }
  // index offset for incrementing z by 1
  m_zoff = static_cast<int>(m_mesh[1].size());
}

auto BFieldMeshZR::memSize() const -> int
{
  // accumulate in size_t; the final conversion to int gives the same value
  // as converting after every addition (modular arithmetic)
  std::size_t size = 0;
  size += sizeof(double) * 6;
  size += sizeof(int) * 1;
  for (int i = 0; i < 2; i++) {
    size += sizeof(double) * m_mesh[i].capacity();
    size += sizeof(int) * m_LUT[i].capacity();
  }
  size += sizeof(BFieldVectorZR) * m_field.capacity();
  return static_cast<int>(size);
}
