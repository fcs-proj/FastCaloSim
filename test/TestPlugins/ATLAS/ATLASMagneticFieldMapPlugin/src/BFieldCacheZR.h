/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

//
// BFieldCacheZR.h
//
// Cache of one bin of the magnetic field map.
// Defined by ranges in z, r, and the (Bz, Br) vectors at the 4 corners of the
// "bin".
//
// Masahiro Morii, Harvard University
//
#ifndef BFIELDCACHEZR_H
#define BFIELDCACHEZR_H

#include <array>
#include <cmath>
#include <cstddef>
#include <iostream>

#include "BFieldVectorZR.h"

class BFieldCacheZR
{
public:
  // default constructor sets unphysical boundaries, so that inside() will fail
  BFieldCacheZR() = default;
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
    m_invz = static_cast<float>(1.0 / (zmax - zmin));
    m_invr = static_cast<float>(1.0 / (rmax - rmin));
  }
  // set the field values at each corner
  void setField(std::size_t i, const BFieldVectorZR& field)
  {
    m_field[0][i] = field[0];
    m_field[1][i] = field[1];
  }
  // set the multiplicative factor for the field vectors
  // test if (z, r) is inside this bin
  auto inside(double z, double r) const -> bool
  {
    return (z >= m_zmin && z <= m_zmax && r >= m_rmin && r <= m_rmax);
  }
  // interpolate the field and return B[3].
  // also compute field derivatives if deriv[9] is given.
  inline void getB(const double* xyz,
                   double r,
                   double* B,
                   double* deriv = nullptr) const;

private:
  double m_zmin {0.0}, m_zmax {-1.0};  // bin range in z
  double m_rmin {0.0}, m_rmax {-1.0};  // bin range in r
  float m_invz {}, m_invr {};  // 1/(bin size) in z, r
  std::array<std::array<float, 4>, 2> m_field {};  // (Bz,Br) at 4 corners
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
  // computed in double and rounded to float; 1.0F - f is identical to
  // float(1.0 - double(f)) (double has more than 2*24+2 mantissa bits)
  auto const fz =
      static_cast<float>((z - m_zmin) * static_cast<double>(m_invz));
  float const gz = 1.0F - fz;
  auto const fr =
      static_cast<float>((r - m_rmin) * static_cast<double>(m_invr));
  float const gr = 1.0F - fr;
  // interpolate field values in z, r
  std::array<float, 2> Bzr {};
  for (std::size_t i = 0; i < 2; i++) {  // z, r components
    const auto& field = m_field[i];
    Bzr[i] = gz * (gr * field[0] + fr * field[1])
        + fz * (gr * field[2] + fr * field[3]);
  }
  // convert (Bz,Br) to (Bx,By,Bz)
  const float invr = r > 0.0 ? static_cast<float>(1.0 / r) : 0.0F;
  auto const c = static_cast<float>(x * static_cast<double>(invr));
  auto const es = static_cast<float>(y * static_cast<double>(invr));
  B[0] = static_cast<double>(Bzr[1] * c);
  B[1] = static_cast<double>(Bzr[1] * es);
  B[2] = static_cast<double>(Bzr[0]);

  // compute field derivatives if requested
  if (deriv != nullptr) {
    std::array<float, 2> dBdz {};
    std::array<float, 2> dBdr {};
    for (std::size_t j = 0; j < 2; j++) {  // Bz, Br components
      const auto& field = m_field[j];
      dBdz[j] =
          m_invz * (gr * (field[2] - field[0]) + fr * (field[3] - field[1]));
      dBdr[j] =
          m_invr * (gz * (field[1] - field[0]) + fz * (field[3] - field[2]));
    }
    // convert to cartesian coordinates
    float const cc = c * c;
    float const cs = c * es;
    float const ss = es * es;
    float const sinvr = es * invr;
    float const cinvr = c * invr;
    deriv[0] =
        static_cast<double>(cc * dBdr[1]) + static_cast<double>(sinvr) * B[1];
    deriv[1] =
        static_cast<double>(cs * dBdr[1]) - static_cast<double>(cinvr) * B[1];
    deriv[2] = static_cast<double>(c * dBdz[1]);
    deriv[3] =
        static_cast<double>(cs * dBdr[1]) - static_cast<double>(sinvr) * B[0];
    deriv[4] =
        static_cast<double>(ss * dBdr[1]) + static_cast<double>(cinvr) * B[0];
    deriv[5] = static_cast<double>(es * dBdz[1]);
    deriv[6] = static_cast<double>(c * dBdr[0]);
    deriv[7] = static_cast<double>(es * dBdr[0]);
    deriv[8] = static_cast<double>(dBdz[0]);
  }
}

#endif
