/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

//
// BFieldMesh.h
//
// Generic 3-d mesh representing a simple field map.
// The field type is templated - it may be short (for the toroid) or double (for
// the solenoid)
//
// Masahiro Morii, Harvard University
//
#ifndef BFIELDMESH_H
#define BFIELDMESH_H

#include <array>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <vector>

#include "BFieldCache.h"
#include "BFieldVector.h"

template<class T>
class BFieldMesh
{
public:
  // constructor
  BFieldMesh()
      : m_scale(1.0)
  {
    ;
  }
  BFieldMesh(double zmin,
             double zmax,
             double rmin,
             double rmax,
             double phimin,
             double phimax,
             double bscale)
      : m_min {zmin, rmin, phimin}
      , m_max {zmax, rmax, phimax}
      , m_scale(bscale)
      , m_nomScale(bscale)
  {
  }
  // set ranges
  void setRange(double zmin,
                double zmax,
                double rmin,
                double rmax,
                double phimin,
                double phimax)
  {
    m_min[0] = zmin;
    m_max[0] = zmax;
    m_min[1] = rmin;
    m_max[1] = rmax;
    m_min[2] = phimin;
    m_max[2] = phimax;
  }
  // set bscale
  void setBscale(double bscale) { m_scale = m_nomScale = bscale; }
  // scale bscale by a factor
  void scaleBscale(double factor) { m_scale = factor * m_nomScale; }
  // allocate space to vectors
  void reserve(int nz, int nr, int nphi, int nfield);
  void reserve(int nz, int nr, int nphi)
  {
    reserve(nz, nr, nphi, nz * nr * nphi);
  }
  // add elements to vectors
  void appendMesh(int i, double mesh)
  {
    m_mesh[static_cast<std::size_t>(i)].push_back(mesh);
  }
  void appendField(const BFieldVector<T>& field) { m_field.push_back(field); }
  // build LUT
  void buildLUT();
  // test if a point is inside this zone
  inline auto inside(double z, double r, double phi) const -> bool;
  // find the bin
  inline void getCache(double z,
                       double r,
                       double phi,
                       BFieldCache& cache) const;
  // get the B field
  void getB(const double* xyz, double* B, double* deriv = nullptr) const;
  // accessors
  auto min(int i) const -> double { return m_min[static_cast<std::size_t>(i)]; }
  auto max(int i) const -> double { return m_max[static_cast<std::size_t>(i)]; }
  auto zmin() const -> double { return m_min[0]; }
  auto zmax() const -> double { return m_max[0]; }
  auto rmin() const -> double { return m_min[1]; }
  auto rmax() const -> double { return m_max[1]; }
  auto phimin() const -> double { return m_min[2]; }
  auto phimax() const -> double { return m_max[2]; }
  auto nmesh(int i) const -> unsigned
  {
    return static_cast<unsigned>(m_mesh[static_cast<std::size_t>(i)].size());
  }
  auto mesh(int i, int j) const -> double
  {
    return m_mesh[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)];
  }
  auto nfield() const -> unsigned
  {
    return static_cast<unsigned>(m_field.size());
  }
  auto field(int i) const -> const BFieldVector<T>&
  {
    return m_field[static_cast<std::size_t>(i)];
  }
  auto bscale() const -> double { return m_scale; }
  // memory used by the mesh (BFieldZone::memSize() adds its conductors)
  auto meshMemSize() const -> int;

private:
  // tentative mesh index in dimension j from the LUT, x = coordinate - min(j)
  auto lutLookup(std::size_t j, double x) const -> std::size_t
  {
    const auto ilut = static_cast<int>(x * m_invUnit[j]);  // index to LUT
    return static_cast<std::size_t>(m_LUT[j][static_cast<std::size_t>(ilut)]);
  }

protected:
  std::array<double, 3> m_min {}, m_max {};
  std::array<std::vector<double>, 3> m_mesh;

private:
  std::vector<BFieldVector<T> > m_field;
  double m_scale;
  double m_nomScale {};  // nominal m_scale from the map
  // look-up table and related variables
  std::array<std::vector<int>, 3> m_LUT;
  std::array<double, 3> m_invUnit {};  // inverse unit size in the LUT
  int m_roff {}, m_zoff {};
};

//
// Implemnetation follows
//

//
// Reserve space in the vectors to avoid unnecessary memory re-allocations.
//
template<class T>
void BFieldMesh<T>::reserve(int nz, int nr, int nphi, int nfield)
{
  m_mesh[0].reserve(static_cast<std::size_t>(nz));
  m_mesh[1].reserve(static_cast<std::size_t>(nr));
  m_mesh[2].reserve(static_cast<std::size_t>(nphi));
  m_field.reserve(static_cast<std::size_t>(nfield));
}

//
// Test if a point (z,r,phi) is inside this mesh region.
//
template<class T>
auto BFieldMesh<T>::inside(double z, double r, double phi) const -> bool
{
  // assume phi is in [-pi,pi].
  // phimax() is in [0,2pi], but phimin() may be < 0 if the range crosses phi =
  // 0. we have to test twice to get all possible cases.
  if (phi < phimin()) {
    phi += 2.0 * M_PI;
  }
  return (phi >= phimin() && phi <= phimax() && z >= zmin() && z <= zmax()
          && r >= rmin() && r <= rmax());
}

//
// Find and return the cache of the bin containing (z,r,phi)
//
template<class T>
void BFieldMesh<T>::getCache(double z,
                             double r,
                             double phi,
                             BFieldCache& cache) const
{
  // make sure phi is inside this zone
  if (phi < phimin()) {
    phi += 2.0 * M_PI;
  }
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
  // phi
  const std::vector<double>& mphi(m_mesh[2]);
  std::size_t iphi =
      lutLookup(2, phi - phimin());  // tentative mesh index from LUT
  if (phi > mphi[iphi + 1]) {
    iphi++;
  }
  // store the bin edges
  cache.setRange(
      mz[iz], mz[iz + 1], mr[ir], mr[ir + 1], mphi[iphi], mphi[iphi + 1]);
  // store the B field at the 8 corners
  const auto roff = static_cast<std::size_t>(m_roff);
  const auto zoff = static_cast<std::size_t>(m_zoff);
  std::size_t const im0 =
      iz * zoff + ir * roff + iphi;  // index of the first corner
  cache.setField(0, m_field[im0]);
  cache.setField(1, m_field[im0 + 1]);
  cache.setField(2, m_field[im0 + roff]);
  cache.setField(3, m_field[im0 + roff + 1]);
  cache.setField(4, m_field[im0 + zoff]);
  cache.setField(5, m_field[im0 + zoff + 1]);
  cache.setField(6, m_field[im0 + zoff + roff]);
  cache.setField(7, m_field[im0 + zoff + roff + 1]);
  // store the B scale
  cache.setBscale(m_scale);
}

//
// Compute the magnetic field (and the derivatives) without caching
//
template<class T>
void BFieldMesh<T>::getB(const double* xyz, double* B, double* deriv) const
{
  // cylindrical coordinates
  double const x = xyz[0];
  double const y = xyz[1];
  double const z = xyz[2];
  double const r = sqrt(x * x + y * y);
  double phi = atan2(y, x);
  if (phi < phimin()) {
    phi += 2.0 * M_PI;
  }
  // is it inside this map?
  if (!inside(z, r, phi)) {  // no
    B[0] = B[1] = B[2] = 0.0;
    if (deriv != nullptr) {
      for (int i = 0; i < 9; i++) {
        deriv[i] = 0.0;
      }
    }
    return;
  }
  // find the bin
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
  // phi
  const std::vector<double>& mphi(m_mesh[2]);
  std::size_t iphi =
      lutLookup(2, phi - phimin());  // tentative mesh index from LUT
  if (phi > mphi[iphi + 1]) {
    iphi++;
  }
  // get the B field at the 8 corners
  const auto roff = static_cast<std::size_t>(m_roff);
  const auto zoff = static_cast<std::size_t>(m_zoff);
  std::size_t const im0 =
      iz * zoff + ir * roff + iphi;  // index of the first corner
  const std::array<BFieldVector<T>, 8> field {m_field[im0],
                                              m_field[im0 + 1],
                                              m_field[im0 + roff],
                                              m_field[im0 + roff + 1],
                                              m_field[im0 + zoff],
                                              m_field[im0 + zoff + 1],
                                              m_field[im0 + zoff + roff],
                                              m_field[im0 + zoff + roff + 1]};
  // fractional position inside this mesh
  double const fz = (z - mz[iz]) / (mz[iz + 1] - mz[iz]);
  double const gz = 1.0 - fz;
  double const fr = (r - mr[ir]) / (mr[ir + 1] - mr[ir]);
  double const gr = 1.0 - fr;
  double const fphi = (phi - mphi[iphi]) / (mphi[iphi + 1] - mphi[iphi]);
  double const gphi = 1.0 - fphi;
  // interpolate field values in z, r, phi
  std::array<double, 3> Bzrphi {};
  for (std::size_t i = 0; i < 3; i++) {  // z, r, phi
    Bzrphi[i] = m_scale
        * (gz
               * (gr * (gphi * field[0][i] + fphi * field[1][i])
                  + fr * (gphi * field[2][i] + fphi * field[3][i]))
           + fz
               * (gr * (gphi * field[4][i] + fphi * field[5][i])
                  + fr * (gphi * field[6][i] + fphi * field[7][i])));
  }
  // convert (Bz,Br,Bphi) to (Bx,By,Bz)
  double const c = (r > 0.0) ? x / r : cos(mphi[iphi]);
  double const es = (r > 0.0) ? y / r : sin(mphi[iphi]);
  B[0] = Bzrphi[1] * c - Bzrphi[2] * es;
  B[1] = Bzrphi[1] * es + Bzrphi[2] * c;
  B[2] = Bzrphi[0];

  // compute field derivatives if requested
  if (deriv != nullptr) {
    double const sz = m_scale / (mz[iz + 1] - mz[iz]);
    double const esr = m_scale / (mr[ir + 1] - mr[ir]);
    double const sphi = m_scale / (mphi[iphi + 1] - mphi[iphi]);
    std::array<double, 3> dBdz {};
    std::array<double, 3> dBdr {};
    std::array<double, 3> dBdphi {};
    for (std::size_t j = 0; j < 3; j++) {  // Bz, Br, Bphi components
      dBdz[j] = sz
          * (gr
                 * (gphi * (field[4][j] - field[0][j])
                    + fphi * (field[5][j] - field[1][j]))
             + fr
                 * (gphi * (field[6][j] - field[2][j])
                    + fphi * (field[7][j] - field[3][j])));
      dBdr[j] = esr
          * (gz
                 * (gphi * (field[2][j] - field[0][j])
                    + fphi * (field[3][j] - field[1][j]))
             + fz
                 * (gphi * (field[6][j] - field[4][j])
                    + fphi * (field[7][j] - field[5][j])));
      dBdphi[j] = sphi
          * (gz
                 * (gr * (field[1][j] - field[0][j])
                    + fr * (field[3][j] - field[2][j]))
             + fz
                 * (gr * (field[5][j] - field[4][j])
                    + fr * (field[7][j] - field[6][j])));
    }
    // convert to cartesian coordinates
    double const cc = c * c;
    double const cs = c * es;
    double const ss = es * es;
    deriv[0] = cc * dBdr[1] - cs * dBdr[2] - cs * dBdphi[1] / r
        + ss * dBdphi[2] / r + es * B[1] / r;
    deriv[1] = cs * dBdr[1] - ss * dBdr[2] + cc * dBdphi[1] / r
        - cs * dBdphi[2] / r - c * B[1] / r;
    deriv[2] = c * dBdz[1] - es * dBdz[2];
    deriv[3] = cs * dBdr[1] + cc * dBdr[2] - ss * dBdphi[1] / r
        - cs * dBdphi[2] / r - es * B[0] / r;
    deriv[5] = es * dBdz[1] + c * dBdz[2];
    deriv[6] = c * dBdr[0] - es * dBdphi[0] / r;
    deriv[7] = es * dBdr[1] + c * dBdphi[0] / r;
    deriv[8] = dBdz[0];
  }
}

//
// Construct the look-up table to accelerate bin-finding.
//
template<class T>
void BFieldMesh<T>::buildLUT()
{
  for (std::size_t j = 0; j < 3; j++) {  // z, r, phi
    // align the m_mesh edges to m_min/m_max
    m_mesh[j].front() = m_min[j];
    m_mesh[j].back() = m_max[j];
    // determine the unit size, q, to be used in the LUTs
    double const width = m_mesh[j].back() - m_mesh[j].front();
    double q(width);
    for (std::size_t i = 0; i < m_mesh[j].size() - 1; i++) {
      q = std::min(q, m_mesh[j][i + 1] - m_mesh[j][i]);
    }
    // find the number of units in the LUT
    int n = static_cast<int>(width / q) + 1;
    q = width / (n + 0.5);
    m_invUnit[j] = 1.0 / q;  // new unit size
    n++;
    int m_number = 0;  // mesh number
    for (int i = 0; i < n; i++) {  // LUT index
      if (i * q + m_mesh[j].front()
          > m_mesh[j][static_cast<std::size_t>(m_number) + 1])
      {
        m_number++;
      }
      m_LUT[j].push_back(m_number);
    }
  }
  // index offset for incrementing r by 1
  m_roff = static_cast<int>(m_mesh[2].size());
  // index offset for incrementing z by 1
  m_zoff =
      static_cast<int>(static_cast<std::size_t>(m_roff) * m_mesh[1].size());
}

template<class T>
auto BFieldMesh<T>::meshMemSize() const -> int
{
  // accumulate in size_t; the final conversion to int gives the same value
  // as converting after every addition (modular arithmetic)
  std::size_t size = 0;
  size += sizeof(double) * 10;
  size += sizeof(int) * 2;
  for (std::size_t i = 0; i < 3; i++) {
    size += sizeof(double) * m_mesh[i].capacity();
    size += sizeof(int) * m_LUT[i].capacity();
  }
  size += sizeof(BFieldVector<T>) * m_field.capacity();
  return static_cast<int>(size);
}

#endif
