/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

//
// BFieldH8Grid.h
//
// A regular x-y-z grid of field data for H8 field map
//
// Masahiro Morii, Harvard University
//
#ifndef BFIELDH8GRID_H
#define BFIELDH8GRID_H

#include <iostream>
#include <vector>

class BFieldH8Grid
{
public:
  // constructor
  BFieldH8Grid();
  // read map data from a text file
  void readMap(std::istream& input);
  // compute magnetic field + derivatives
  void getB(const double* xyz, double* B, double* deriv = 0) const;
  // true if the grid has been defined
  bool defined() const { return (m_n[0] > 0); }
  // true if xyz[3] is inside this grid
  bool inside(const double* xyz) const
  {
    return (xyz[0] >= m_min[0] && xyz[0] <= m_max[0] && xyz[1] >= m_min[1]
            && xyz[1] <= m_max[1] && xyz[2] >= m_min[2] && xyz[2] <= m_max[2]);
  }
  // set x-y-z offset of the grid coordinates
  void setOffset(const double* dxyz);
  // get copies of size vectors
  void getBounds(double* out_min, double* out_max, double* out_d) const
  {
    for (unsigned i = 0; i < 3; i++) {
      out_min[i] = m_min[i];
      out_max[i] = m_max[i];
      out_d[i] = m_d[i];
    }
    return;
  }

private:
  int m_n[3];  // number of grid points
  double m_min[3], m_max[3];  // range in x,y,z (mm)
  double m_d[3];  // offset in x,y,z (mm)
  std::vector<double> m_B[3];  // Bz,By,Bz (kT)
};

#endif
