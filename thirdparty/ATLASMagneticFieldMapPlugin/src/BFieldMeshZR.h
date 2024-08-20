/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

//
// BFieldMeshZR.h
//
// A 2-dim z-r mesh inside the solenoid field map
//
// Masahiro Morii, Harvard University
//
#ifndef BFIELDMESHZR_H
#define BFIELDMESHZR_H

#include <cmath>
#include <vector>

#include "BFieldCacheZR.h"
#include "BFieldVectorZR.h"

class BFieldMeshZR
{
public:
  // constructor
  BFieldMeshZR(double zmin, double zmax, double rmin, double rmax)
  {
    m_min[0] = zmin;
    m_max[0] = zmax;
    m_min[1] = rmin;
    m_max[1] = rmax;
  }
  // allocate space to vectors
  void reserve(int nz, int nr)
  {
    m_mesh[0].reserve(nz);
    m_mesh[1].reserve(nr);
    m_field.reserve(nz * nr);
  }
  // add elements to vectors
  void appendMesh(int i, double mesh) { m_mesh[i].push_back(mesh); }
  void appendField(const BFieldVectorZR& field) { m_field.push_back(field); }
  // build LUT
  void buildLUT();
  // test if a point is inside this zone
  inline bool inside(double z, double r) const
  {
    return (z >= zmin() && z <= zmax() && r >= rmin() && r <= rmax());
  }
  // find the bin
  inline void getCache(double z, double r, BFieldCacheZR& cache) const;
  // accessors
  double min(int i) const { return m_min[i]; }
  double max(int i) const { return m_max[i]; }
  double zmin() const { return m_min[0]; }
  double zmax() const { return m_max[0]; }
  double rmin() const { return m_min[1]; }
  double rmax() const { return m_max[1]; }
  unsigned nmesh(int i) const { return m_mesh[i].size(); }
  double mesh(int i, int j) const { return m_mesh[i][j]; }
  unsigned nfield() const { return m_field.size(); }
  const BFieldVectorZR& field(int i) const { return m_field[i]; }
  int memSize() const;

private:
  double m_min[2], m_max[2];
  std::vector<double> m_mesh[2];
  std::vector<BFieldVectorZR> m_field;
  // look-up table and related variables
  std::vector<int> m_LUT[2];
  double m_invUnit[2];  // inverse unit size in the LUT
  int m_zoff;
};

//
// Find and return the cache of the bin containing (z,r)
//
inline void BFieldMeshZR::getCache(double z,
                                   double r,
                                   BFieldCacheZR& cache) const
{
  // find the mesh, and relative location in the mesh
  // z
  const std::vector<double>& mz(m_mesh[0]);
  int iz = int((z - zmin()) * m_invUnit[0]);  // index to LUT
  iz = m_LUT[0][iz];  // tentative mesh index from LUT
  if (z > mz[iz + 1])
    iz++;
  // r
  const std::vector<double>& mr(m_mesh[1]);
  int ir = int((r - rmin()) * m_invUnit[1]);  // index to LUT
  ir = m_LUT[1][ir];  // tentative mesh index from LUT
  if (r > mr[ir + 1])
    ir++;
  // store the bin edges
  cache.setRange(mz[iz], mz[iz + 1], mr[ir], mr[ir + 1]);
  // store the B field at the 8 corners
  int im0 = iz * m_zoff + ir;  // index of the first corner
  cache.setField(0, m_field[im0]);
  cache.setField(1, m_field[im0 + 1]);
  cache.setField(2, m_field[im0 + m_zoff]);
  cache.setField(3, m_field[im0 + m_zoff + 1]);
  return;
}

#endif
