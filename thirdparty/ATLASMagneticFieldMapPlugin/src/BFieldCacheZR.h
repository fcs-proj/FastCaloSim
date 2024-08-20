/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

//
// BFieldCacheZR.h
//
// Cashe of one bin of the magnetic field map.
// Defined by ranges in z, r, and the (Bz, Br) vectors at the 4 corners of the
// "bin".
//
// Masahiro Morii, Harvard University
//
#ifndef BFIELDCACHEZR_H
#define BFIELDCACHEZR_H

#include <cmath>
#include <iostream>

#include "BFieldVectorZR.h"

class BFieldCacheZR
{
public:
  // default constructor sets unphysical boundaries, so that inside() will fail
  BFieldCacheZR()
      : m_zmin(0.0)
      , m_zmax(-1.0)
      , m_rmin(0.0)
      , m_rmax(-1.0)
  {
    ;
  }
  // invalidate this cache, so that inside() will fail
  void invalidate()
  {
    m_rmin = 0.0;
    m_rmax = -1.0;
  }
  // set the z, r range that defines the bin
  void setRange(double zmin, double zmax, double rmin, double rmax)
  {
    m_zmin = zmin;
    m_zmax = zmax;
    m_rmin = rmin;
    m_rmax = rmax;
    m_invz = 1.0 / (zmax - zmin);
    m_invr = 1.0 / (rmax - rmin);
  }
  // set the field values at each corner
  void setField(int i, const BFieldVectorZR& field)
  {
    m_field[0][i] = field[0];
    m_field[1][i] = field[1];
  }
  // set the multiplicative factor for the field vectors
  // test if (z, r) is inside this bin
  bool inside(double z, double r) const
  {
    return (z >= m_zmin && z <= m_zmax && r >= m_rmin && r <= m_rmax);
  }
  // interpolate the field and return B[3].
  // also compute field derivatives if deriv[9] is given.
  inline void getB(const double* xyz,
                   double r,
                   double* B,
                   double* deriv = 0) const;

private:
  double m_zmin, m_zmax;  // bin range in z
  double m_rmin, m_rmax;  // bin range in r
  float m_invz, m_invr;  // 1/(bin size) in z, r
  float m_field[2][4];  // (Bz,Br) at 4 corners of the bin
};

// inline methods

void BFieldCacheZR::getB(const double* xyz,
                         double r,
                         double* B,
                         double* deriv) const
{
  const double& x(xyz[0]);
  const double& y(xyz[1]);
  const double& z(xyz[2]);
  // fractional position inside this bin
  float fz = (z - m_zmin) * m_invz;
  float gz = 1.0 - fz;
  float fr = (r - m_rmin) * m_invr;
  float gr = 1.0 - fr;
  // interpolate field values in z, r
  float Bzr[2];
  for (int i = 0; i < 2; i++) {  // z, r components
    const float* field = m_field[i];
    Bzr[i] = gz * (gr * field[0] + fr * field[1])
        + fz * (gr * field[2] + fr * field[3]);
  }
  // convert (Bz,Br) to (Bx,By,Bz)
  float invr;
  if (r > 0.0) {
    invr = 1.0 / r;
  } else {
    invr = 0.0;
  }
  float c(x * invr);
  float es(y * invr);
  B[0] = Bzr[1] * c;
  B[1] = Bzr[1] * es;
  B[2] = Bzr[0];

  // compute field derivatives if requested
  if (deriv) {
    float dBdz[2], dBdr[2];
    for (int j = 0; j < 2; j++) {  // Bz, Br components
      const float* field = m_field[j];
      dBdz[j] =
          m_invz * (gr * (field[2] - field[0]) + fr * (field[3] - field[1]));
      dBdr[j] =
          m_invr * (gz * (field[1] - field[0]) + fz * (field[3] - field[2]));
    }
    // convert to cartesian coordinates
    float cc = c * c;
    float cs = c * es;
    float ss = es * es;
    float sinvr = es * invr;
    float cinvr = c * invr;
    deriv[0] = cc * dBdr[1] + sinvr * B[1];
    deriv[1] = cs * dBdr[1] - cinvr * B[1];
    deriv[2] = c * dBdz[1];
    deriv[3] = cs * dBdr[1] - sinvr * B[0];
    deriv[4] = ss * dBdr[1] + cinvr * B[0];
    deriv[5] = es * dBdz[1];
    deriv[6] = c * dBdr[0];
    deriv[7] = es * dBdr[0];
    deriv[8] = dBdz[0];
  }
}

#endif
