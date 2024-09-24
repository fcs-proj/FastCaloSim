/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

//
// BFieldH8Grid.cxx
//
// Masahiro Morii, Harvard University
//
#include "BFieldH8Grid.h"

BFieldH8Grid::BFieldH8Grid()
{
  for (int i = 0; i < 3; i++) {
    m_n[i] = 0;
    m_min[i] = 0.0;
    m_max[i] = 0.0;
    m_d[i] = 0.0;
  }
}

void BFieldH8Grid::readMap(std::istream& input)
{
  // read the magnet header line
  char name[256];
  double cor;
  m_n[0] = 0;
  input >> name;
  for (int i = 0; i < 3; i++)
    input >> m_n[i];
  for (int i = 0; i < 3; i++)
    input >> m_min[i] >> m_max[i];
  for (int i = 0; i < 3; i++)
    input >> m_d[i];
  input >> cor;
  // anything read?
  if (m_n[0] == 0)
    return;
  // convert unit to mm
  const double meter(1000.);  // m in mm
  for (int i = 0; i < 3; i++) {
    m_min[i] *= meter;
    m_max[i] *= meter;
    m_d[i] *= meter;
  }
  // prepare space for the field data
  int nxyz = m_n[0] * m_n[1] * m_n[2];
  for (int i = 0; i < 3; i++)
    m_B[i].resize(nxyz);
  // read the field data
  double xyz[3];
  double B[3];
  for (int k = 0; k < nxyz; k++) {
    input >> xyz[0] >> xyz[1] >> xyz[2] >> B[0] >> B[1] >> B[2];
    // adjust unit and scale
    const double cm(10.);  // cm in mm
    const double tesla(0.001);  // T in kT
    for (int i = 0; i < 3; i++) {
      xyz[i] *= cm;
      B[i] *= -cor * tesla;
    }
    // find the grid number
    int j[3];
    for (int i = 0; i < 3; i++) {
      j[i] =
          int((xyz[i] - m_min[i]) / (m_max[i] - m_min[i]) * (m_n[i] - 1) + 0.5);
    }
    // store
    int ixyz = j[0] + m_n[0] * (j[1] + m_n[1] * j[2]);
    for (int i = 0; i < 3; i++) {
      m_B[i][ixyz] = B[i];
    }
  }
  // apply offset to the range
  for (int i = 0; i < 3; i++) {
    m_min[i] += m_d[i];
    m_max[i] += m_d[i];
  }
}

void BFieldH8Grid::getB(const double* xyz, double* B, double* deriv) const
{
  if (!defined() || !inside(xyz)) {
    B[0] = B[1] = B[2] = 0.0;
    if (deriv != 0) {
      for (int j = 0; j < 9; j++)
        deriv[j] = 0.0;
    }
    return;
  }
  // find the grid index
  double invunit[3];
  int j[3];
  double f[3];
  double g[3];
  for (int i = 0; i < 3; i++) {
    invunit[i] = (m_n[i] - 1) / (m_max[i] - m_min[i]);
    double a = (xyz[i] - m_min[i]) * invunit[i];
    j[i] = int(a);
    f[i] = a - j[i];
    g[i] = 1.0 - f[i];
  }
  int ixyz[8];
  ixyz[0] = j[0] + m_n[0] * (j[1] + m_n[1] * j[2]);
  ixyz[1] = ixyz[0] + 1;
  ixyz[2] = ixyz[0] + m_n[0];
  ixyz[3] = ixyz[2] + 1;
  ixyz[4] = ixyz[0] + m_n[0] * m_n[1];
  ixyz[5] = ixyz[4] + 1;
  ixyz[6] = ixyz[4] + m_n[0];
  ixyz[7] = ixyz[6] + 1;
  // interpolate field values
  for (int i = 0; i < 3; i++) {  // Bx, By, Bz
    B[i] = g[2]
            * (g[1] * (g[0] * m_B[i][ixyz[0]] + f[0] * m_B[i][ixyz[1]])
               + f[1] * (g[0] * m_B[i][ixyz[2]] + f[0] * m_B[i][ixyz[3]]))
        + f[2]
            * (g[1] * (g[0] * m_B[i][ixyz[4]] + f[0] * m_B[i][ixyz[5]])
               + f[1] * (g[0] * m_B[i][ixyz[6]] + f[0] * m_B[i][ixyz[7]]));
  }
  // derivatives
  if (deriv != 0) {
    for (int i = 0; i < 3; i++) {  // Bx, By, Bz
      deriv[i * 3] = (g[2]
                          * (g[1] * (m_B[i][ixyz[1]] - m_B[i][ixyz[0]])
                             + f[1] * (m_B[i][ixyz[3]] - m_B[i][ixyz[2]]))
                      + f[2]
                          * (g[1] * (m_B[i][ixyz[5]] - m_B[i][ixyz[4]])
                             + f[1] * (m_B[i][ixyz[7]] - m_B[i][ixyz[6]])))
          * invunit[0];
      deriv[i * 3 + 1] = (g[2]
                              * (g[0] * (m_B[i][ixyz[2]] - m_B[i][ixyz[0]])
                                 + f[0] * (m_B[i][ixyz[3]] - m_B[i][ixyz[1]]))
                          + f[2]
                              * (g[0] * (m_B[i][ixyz[6]] - m_B[i][ixyz[4]])
                                 + f[0] * (m_B[i][ixyz[7]] - m_B[i][ixyz[5]])))
          * invunit[1];
      deriv[i * 3 + 2] = (g[1]
                              * (g[0] * (m_B[i][ixyz[4]] - m_B[i][ixyz[0]])
                                 + f[0] * (m_B[i][ixyz[5]] - m_B[i][ixyz[1]]))
                          + f[1]
                              * (g[0] * (m_B[i][ixyz[6]] - m_B[i][ixyz[2]])
                                 + f[0] * (m_B[i][ixyz[7]] - m_B[i][ixyz[3]])))
          * invunit[2];
    }
  }
}

void BFieldH8Grid::setOffset(const double* dxyz)
{
  // update the range by the change in offset
  for (int i = 0; i < 3; i++) {
    m_min[i] += dxyz[i] - m_d[i];
    m_max[i] += dxyz[i] - m_d[i];
    m_d[i] = dxyz[i];
  }
}
