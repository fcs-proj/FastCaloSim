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

#include <cmath>
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
      : m_scale(bscale)
      , m_nomScale(bscale)
  {
    m_min[0] = zmin;
    m_max[0] = zmax;
    m_min[1] = rmin;
    m_max[1] = rmax;
    m_min[2] = phimin;
    m_max[2] = phimax;
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
  void appendMesh(int i, double mesh) { m_mesh[i].push_back(mesh); }
  void appendField(const BFieldVector<T>& field) { m_field.push_back(field); }
  // build LUT
  void buildLUT();
  // test if a point is inside this zone
  inline bool inside(double z, double r, double phi) const;
  // find the bin
  inline void getCache(double z,
                       double r,
                       double phi,
                       BFieldCache& cache) const;
  // get the B field
  void getB(const double* xyz, double* B, double* deriv = 0) const;
  // accessors
  double min(int i) const { return m_min[i]; }
  double max(int i) const { return m_max[i]; }
  double zmin() const { return m_min[0]; }
  double zmax() const { return m_max[0]; }
  double rmin() const { return m_min[1]; }
  double rmax() const { return m_max[1]; }
  double phimin() const { return m_min[2]; }
  double phimax() const { return m_max[2]; }
  unsigned nmesh(int i) const { return m_mesh[i].size(); }
  double mesh(int i, int j) const { return m_mesh[i][j]; }
  unsigned nfield() const { return m_field.size(); }
  const BFieldVector<T>& field(int i) const { return m_field[i]; }
  double bscale() const { return m_scale; }
  int memSize() const;

protected:
  double m_min[3], m_max[3];
  std::vector<double> m_mesh[3];

private:
  std::vector<BFieldVector<T> > m_field;
  double m_scale;
  double m_nomScale;  // nominal m_scale from the map
  // look-up table and related variables
  std::vector<int> m_LUT[3];
  double m_invUnit[3];  // inverse unit size in the LUT
  int m_roff, m_zoff;
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
  m_mesh[0].reserve(nz);
  m_mesh[1].reserve(nr);
  m_mesh[2].reserve(nphi);
  m_field.reserve(nfield);
}

//
// Test if a point (z,r,phi) is inside this mesh region.
//
template<class T>
bool BFieldMesh<T>::inside(double z, double r, double phi) const
{
  // assume phi is in [-pi,pi].
  // phimax() is in [0,2pi], but phimin() may be < 0 if the range crosses phi =
  // 0. we have to test twice to get all possible cases.
  if (phi < phimin())
    phi += 2.0 * M_PI;
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
  if (phi < phimin())
    phi += 2.0 * M_PI;
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
  // phi
  const std::vector<double>& mphi(m_mesh[2]);
  int iphi = int((phi - phimin()) * m_invUnit[2]);  // index to LUT
  iphi = m_LUT[2][iphi];  // tentative mesh index from LUT
  if (phi > mphi[iphi + 1])
    iphi++;
  // store the bin edges
  cache.setRange(
      mz[iz], mz[iz + 1], mr[ir], mr[ir + 1], mphi[iphi], mphi[iphi + 1]);
  // store the B field at the 8 corners
  int im0 = iz * m_zoff + ir * m_roff + iphi;  // index of the first corner
  cache.setField(0, m_field[im0]);
  cache.setField(1, m_field[im0 + 1]);
  cache.setField(2, m_field[im0 + m_roff]);
  cache.setField(3, m_field[im0 + m_roff + 1]);
  cache.setField(4, m_field[im0 + m_zoff]);
  cache.setField(5, m_field[im0 + m_zoff + 1]);
  cache.setField(6, m_field[im0 + m_zoff + m_roff]);
  cache.setField(7, m_field[im0 + m_zoff + m_roff + 1]);
  // store the B scale
  cache.setBscale(m_scale);
  return;
}

//
// Compute the magnetic field (and the derivatives) without caching
//
template<class T>
void BFieldMesh<T>::getB(const double* xyz, double* B, double* deriv) const
{
  // cylindrical coordinates
  double x = xyz[0];
  double y = xyz[1];
  double z = xyz[2];
  double r = sqrt(x * x + y * y);
  double phi = atan2(y, x);
  if (phi < phimin())
    phi += 2.0 * M_PI;
  // is it inside this map?
  if (!inside(z, r, phi)) {  // no
    B[0] = B[1] = B[2] = 0.0;
    if (deriv)
      for (int i = 0; i < 9; i++)
        deriv[i] = 0.0;
    return;
  }
  // find the bin
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
  // phi
  const std::vector<double>& mphi(m_mesh[2]);
  int iphi = int((phi - phimin()) * m_invUnit[2]);  // index to LUT
  iphi = m_LUT[2][iphi];  // tentative mesh index from LUT
  if (phi > mphi[iphi + 1])
    iphi++;
  // get the B field at the 8 corners
  int im0 = iz * m_zoff + ir * m_roff + iphi;  // index of the first corner
  BFieldVector<T> field[8];
  field[0] = m_field[im0];
  field[1] = m_field[im0 + 1];
  field[2] = m_field[im0 + m_roff];
  field[3] = m_field[im0 + m_roff + 1];
  field[4] = m_field[im0 + m_zoff];
  field[5] = m_field[im0 + m_zoff + 1];
  field[6] = m_field[im0 + m_zoff + m_roff];
  field[7] = m_field[im0 + m_zoff + m_roff + 1];
  // fractional position inside this mesh
  double fz = (z - mz[iz]) / (mz[iz + 1] - mz[iz]);
  double gz = 1.0 - fz;
  double fr = (r - mr[ir]) / (mr[ir + 1] - mr[ir]);
  double gr = 1.0 - fr;
  double fphi = (phi - mphi[iphi]) / (mphi[iphi + 1] - mphi[iphi]);
  double gphi = 1.0 - fphi;
  // interpolate field values in z, r, phi
  double Bzrphi[3];
  for (int i = 0; i < 3; i++) {  // z, r, phi
    Bzrphi[i] = m_scale
        * (gz
               * (gr * (gphi * field[0][i] + fphi * field[1][i])
                  + fr * (gphi * field[2][i] + fphi * field[3][i]))
           + fz
               * (gr * (gphi * field[4][i] + fphi * field[5][i])
                  + fr * (gphi * field[6][i] + fphi * field[7][i])));
  }
  // convert (Bz,Br,Bphi) to (Bx,By,Bz)
  double c = (r > 0.0) ? x / r : cos(mphi[iphi]);
  double es = (r > 0.0) ? y / r : sin(mphi[iphi]);
  B[0] = Bzrphi[1] * c - Bzrphi[2] * es;
  B[1] = Bzrphi[1] * es + Bzrphi[2] * c;
  B[2] = Bzrphi[0];

  // compute field derivatives if requested
  if (deriv) {
    double sz = m_scale / (mz[iz + 1] - mz[iz]);
    double esr = m_scale / (mr[ir + 1] - mr[ir]);
    double sphi = m_scale / (mphi[iphi + 1] - mphi[iphi]);
    double dBdz[3], dBdr[3], dBdphi[3];
    for (int j = 0; j < 3; j++) {  // Bz, Br, Bphi components
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
    double cc = c * c;
    double cs = c * es;
    double ss = es * es;
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
  return;
}

//
// Construct the look-up table to accelerate bin-finding.
//
template<class T>
void BFieldMesh<T>::buildLUT()
{
  for (int j = 0; j < 3; j++) {  // z, r, phi
    // align the m_mesh edges to m_min/m_max
    m_mesh[j].front() = m_min[j];
    m_mesh[j].back() = m_max[j];
    // determine the unit size, q, to be used in the LUTs
    double width = m_mesh[j].back() - m_mesh[j].front();
    double q(width);
    for (unsigned i = 0; i < m_mesh[j].size() - 1; i++) {
      q = std::min(q, m_mesh[j][i + 1] - m_mesh[j][i]);
    }
    // find the number of units in the LUT
    int n = int(width / q) + 1;
    q = width / (n + 0.5);
    m_invUnit[j] = 1.0 / q;  // new unit size
    n++;
    int m_number = 0;  // mesh number
    for (int i = 0; i < n; i++) {  // LUT index
      if (i * q + m_mesh[j].front() > m_mesh[j][m_number + 1])
        m_number++;
      m_LUT[j].push_back(m_number);
    }
  }
  m_roff = m_mesh[2].size();  // index offset for incrementing r by 1
  m_zoff = m_roff * m_mesh[1].size();  // index offset for incrementing z by 1
}

template<class T>
int BFieldMesh<T>::memSize() const
{
  int size = 0;
  size += sizeof(double) * 10;
  size += sizeof(int) * 2;
  for (int i = 0; i < 3; i++) {
    size += sizeof(double) * m_mesh[i].capacity();
    size += sizeof(int) * m_LUT[i].capacity();
  }
  size += sizeof(BFieldVector<T>) * m_field.capacity();
  return size;
}

#endif
