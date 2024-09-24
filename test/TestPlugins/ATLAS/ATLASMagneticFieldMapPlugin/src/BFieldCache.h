/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

//
// BFieldCache.h
//
// Cache of one bin of the magnetic field map.
// Defined by ranges in z, r, phi, and the B vectors at the 8 corners of the
// "bin".
//
// Masahiro Morii, Harvard University
//
#ifndef BFIELDCACHE_H
#define BFIELDCACHE_H

#include <cmath>
#include <iostream>

#include "BFieldVector.h"

class BFieldCache
{
public:
  // default constructor sets unphysical boundaries, so that inside() will fail
  BFieldCache()
      : m_zmin(0.0)
      , m_zmax(0.0)
      , m_rmin(0.0)
      , m_rmax(0.0)
      , m_phimin(0.0)
      , m_phimax(-1.0)
  {
    ;
  }
  // make this cache invalid, so that inside() will fail
  void invalidate()
  {
    m_phimin = 0.0;
    m_phimax = -1.0;
  }
  // set the z, r, phi range that defines the bin
  void setRange(double zmin,
                double zmax,
                double rmin,
                double rmax,
                double phimin,
                double phimax)
  {
    m_zmin = zmin;
    m_zmax = zmax;
    m_rmin = rmin;
    m_rmax = rmax;
    m_phimin = phimin;
    m_phimax = phimax;
    m_invz = 1.0 / (zmax - zmin);
    m_invr = 1.0 / (rmax - rmin);
    m_invphi = 1.0 / (phimax - phimin);
  }
  // set the field values at each corner
  void setField(int i, const BFieldVector<double>& field)
  {
    for (int j = 0; j < 3; j++)
      m_field[j][i] = field[j];
  }
  void setField(int i, const BFieldVector<short>& field)
  {
    for (int j = 0; j < 3; j++)
      m_field[j][i] = field[j];
  }  // set the multiplicative factor for the field vectors
  void setBscale(double bscale) { m_scale = bscale; }
  // test if (z, r, phi) is inside this bin
  bool inside(double z, double r, double phi) const
  {
    if (phi < m_phimin)
      phi += 2.0 * M_PI;
    return (phi >= m_phimin && phi <= m_phimax && z >= m_zmin && z <= m_zmax
            && r >= m_rmin && r <= m_rmax);
  }
  // interpolate the field and return B[3].
  // also compute field derivatives if deriv[9] is given.
  inline void getB(const double* xyz,
                   double r,
                   double phi,
                   double* B,
                   double* deriv = 0) const;

private:
  double m_zmin, m_zmax;  // bin range in z
  double m_rmin, m_rmax;  // bin range in r
  double m_phimin, m_phimax;  // bin range in phi
  float m_invz, m_invr, m_invphi;  // 1/(bin size) in z, r, phi
  float m_field[3][8];  // (Bz,Br,Bphi) at 8 corners of the bin
  float m_scale;  // unit of m_field in kT
};

// inline methods

void BFieldCache::getB(
    const double* xyz, double r, double phi, double* B, double* deriv) const
{
  const double& x(xyz[0]);
  const double& y(xyz[1]);
  const double& z(xyz[2]);
  // make sure phi is inside [m_phimin,m_phimax]
  if (phi < m_phimin)
    phi += 2 * M_PI;
  // fractional position inside this bin
  float fz = (z - m_zmin) * m_invz;
  float gz = 1.0 - fz;
  float fr = (r - m_rmin) * m_invr;
  float gr = 1.0 - fr;
  float fphi = (phi - m_phimin) * m_invphi;
  float gphi = 1.0 - fphi;
  // interpolate field values in z, r, phi
  float Bzrphi[3];
  for (int i = 0; i < 3; i++) {  // z, r, phi components
    const float* field = m_field[i];
    Bzrphi[i] = m_scale
        * (gz
               * (gr * (gphi * field[0] + fphi * field[1])
                  + fr * (gphi * field[2] + fphi * field[3]))
           + fz
               * (gr * (gphi * field[4] + fphi * field[5])
                  + fr * (gphi * field[6] + fphi * field[7])));
  }
  // convert (Bz,Br,Bphi) to (Bx,By,Bz)
  float invr, c, es;
  if (r > 0.0) {
    invr = 1.0 / r;
    c = x * invr;
    es = y * invr;
  } else {
    invr = 0.0;
    c = cos(m_phimin);
    es = sin(m_phimin);
  }
  B[0] = Bzrphi[1] * c - Bzrphi[2] * es;
  B[1] = Bzrphi[1] * es + Bzrphi[2] * c;
  B[2] = Bzrphi[0];

  // compute field derivatives if requested
  if (deriv) {
    float sz = m_scale * m_invz;
    float esr = m_scale * m_invr;
    float sphi = m_scale * m_invphi;
    float dBdz[3], dBdr[3], dBdphi[3];
    for (int j = 0; j < 3; j++) {  // Bz, Br, Bphi components
      const float* field = m_field[j];
      dBdz[j] = sz
          * (gr * (gphi * (field[4] - field[0]) + fphi * (field[5] - field[1]))
             + fr
                 * (gphi * (field[6] - field[2])
                    + fphi * (field[7] - field[3])));
      dBdr[j] = esr
          * (gz * (gphi * (field[2] - field[0]) + fphi * (field[3] - field[1]))
             + fz
                 * (gphi * (field[6] - field[4])
                    + fphi * (field[7] - field[5])));
      dBdphi[j] = sphi
          * (gz * (gr * (field[1] - field[0]) + fr * (field[3] - field[2]))
             + fz * (gr * (field[5] - field[4]) + fr * (field[7] - field[6])));
    }
    // convert to cartesian coordinates
    float cc = c * c;
    float cs = c * es;
    float ss = es * es;
    float ccinvr = cc * invr;
    float csinvr = cs * invr;
    float ssinvr = ss * invr;
    float sinvr = es * invr;
    float cinvr = c * invr;
    deriv[0] = cc * dBdr[1] - cs * dBdr[2] - csinvr * dBdphi[1]
        + ssinvr * dBdphi[2] + sinvr * B[1];
    deriv[1] = cs * dBdr[1] - ss * dBdr[2] + ccinvr * dBdphi[1]
        - csinvr * dBdphi[2] - cinvr * B[1];
    deriv[2] = c * dBdz[1] - es * dBdz[2];
    deriv[3] = cs * dBdr[1] + cc * dBdr[2] - ssinvr * dBdphi[1]
        - csinvr * dBdphi[2] - sinvr * B[0];
    deriv[4] = ss * dBdr[1] + cs * dBdr[2] + csinvr * dBdphi[1]
        + ccinvr * dBdphi[2] + cinvr * B[0];
    deriv[5] = es * dBdz[1] + c * dBdz[2];
    deriv[6] = c * dBdr[0] - sinvr * dBdphi[0];
    deriv[7] = es * dBdr[0] + cinvr * dBdphi[0];
    deriv[8] = dBdz[0];
  }
}

#endif
