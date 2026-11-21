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
#include <cstddef>
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
    m_mesh[0].reserve(static_cast<std::size_t>(nz));
    m_mesh[1].reserve(static_cast<std::size_t>(nr));
    m_field.reserve(static_cast<std::size_t>(nz * nr));
  }
  // add elements to vectors
  void appendMesh(int i, double mesh) { m_mesh[i].push_back(mesh); }
  void appendField(const BFieldVectorZR& field) { m_field.push_back(field); }
  // build LUT
  void buildLUT();
  // test if a point is inside this zone
  auto inside(double z, double r) const -> bool
  {
    return (z >= zmin() && z <= zmax() && r >= rmin() && r <= rmax());
  }
  // find the bin
  inline void getCache(double z, double r, BFieldCacheZR& cache) const;
  // accessors
  auto min(int i) const -> double { return m_min[i]; }
  auto max(int i) const -> double { return m_max[i]; }
  auto zmin() const -> double { return m_min[0]; }
  auto zmax() const -> double { return m_max[0]; }
  auto rmin() const -> double { return m_min[1]; }
  auto rmax() const -> double { return m_max[1]; }
  auto nmesh(int i) const -> unsigned
  {
    return static_cast<unsigned>(m_mesh[i].size());
  }
  auto mesh(int i, int j) const -> double
  {
    return m_mesh[i][static_cast<std::size_t>(j)];
  }
  auto nfield() const -> unsigned
  {
    return static_cast<unsigned>(m_field.size());
  }
  auto field(int i) const -> const BFieldVectorZR&
  {
    return m_field[static_cast<std::size_t>(i)];
  }
  auto memSize() const -> int;

private:
  // tentative mesh index in dimension j from the LUT, x = coordinate - min(j)
  auto lutLookup(int j, double x) const -> std::size_t
  {
    const auto ilut = static_cast<int>(x * m_invUnit[j]);  // index to LUT
    return static_cast<std::size_t>(m_LUT[j][static_cast<std::size_t>(ilut)]);
  }

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
  std::size_t iz = lutLookup(0, z - zmin());  // tentative mesh index from LUT
  if (z > mz[iz + 1]) {
    iz++;
  }
  // r
  const std::vector<double>& mr(m_mesh[1]);
  std::size_t ir = lutLookup(1, r - rmin());  // tentative mesh index from LUT
  if (r > mr[ir + 1]) {
    ir++;
  }
  // store the bin edges
  cache.setRange(mz[iz], mz[iz + 1], mr[ir], mr[ir + 1]);
  // store the B field at the 8 corners
  const auto zoff = static_cast<std::size_t>(m_zoff);
  std::size_t const im0 = iz * zoff + ir;  // index of the first corner
  cache.setField(0, m_field[im0]);
  cache.setField(1, m_field[im0 + 1]);
  cache.setField(2, m_field[im0 + zoff]);
  cache.setField(3, m_field[im0 + zoff + 1]);
}

#endif
