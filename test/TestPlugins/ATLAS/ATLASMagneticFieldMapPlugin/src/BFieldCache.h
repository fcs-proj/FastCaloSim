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

#include <array>
#include <cmath>
#include <cstddef>
#include <iostream>

#include "BFieldVector.h"

class BFieldCache
{
public:
  // default constructor sets unphysical boundaries, so that inside() will fail
  BFieldCache() = default;
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
    m_invz = static_cast<float>(1.0 / (zmax - zmin));
    m_invr = static_cast<float>(1.0 / (rmax - rmin));
    m_invphi = static_cast<float>(1.0 / (phimax - phimin));
  }
  // set the field values at each corner
  void setField(std::size_t i, const BFieldVector<double>& field)
  {
    for (std::size_t j = 0; j < 3; j++) {
      m_field[j][i] = static_cast<float>(field[j]);
    }
  }
  void setField(std::size_t i, const BFieldVector<short>& field)
  {
    for (std::size_t j = 0; j < 3; j++) {
      m_field[j][i] = field[j];
    }
  }  // set the multiplicative factor for the field vectors
  void setBscale(double bscale) { m_scale = static_cast<float>(bscale); }
  // test if (z, r, phi) is inside this bin
  auto inside(double z, double r, double phi) const -> bool
  {
    if (phi < m_phimin) {
      phi += 2.0 * M_PI;
    }
    return (phi >= m_phimin && phi <= m_phimax && z >= m_zmin && z <= m_zmax
            && r >= m_rmin && r <= m_rmax);
  }
  // interpolate the field and return B[3].
  // also compute field derivatives if deriv[9] is given.
  inline void getB(const double* xyz,
                   double r,
                   double phi,
                   double* B,
                   double* deriv = nullptr) const;

private:
  double m_zmin {0.0}, m_zmax {0.0};  // bin range in z
  double m_rmin {0.0}, m_rmax {0.0};  // bin range in r
  double m_phimin {0.0}, m_phimax {-1.0};  // bin range in phi
  float m_invz {}, m_invr {}, m_invphi {};  // 1/(bin size) in z, r, phi
  // (Bz,Br,Bphi) at 8 corners of the bin
  std::array<std::array<float, 8>, 3> m_field {};
  float m_scale {};  // unit of m_field in kT
};

// inline methods

void BFieldCache::getB(
    const double* xyz, double r, double phi, double* B, double* deriv) const
{
  const double& x(xyz[0]);
  const double& y(xyz[1]);
  const double& z(xyz[2]);
  // make sure phi is inside [m_phimin,m_phimax]
  if (phi < m_phimin) {
    phi += 2 * M_PI;
  }
  // fractional position inside this bin
  // computed in double and rounded to float; 1.0F - f is identical to
  // float(1.0 - double(f)) (double has more than 2*24+2 mantissa bits)
  auto const fz =
      static_cast<float>((z - m_zmin) * static_cast<double>(m_invz));
  float const gz = 1.0F - fz;
  auto const fr =
      static_cast<float>((r - m_rmin) * static_cast<double>(m_invr));
  float const gr = 1.0F - fr;
  auto const fphi =
      static_cast<float>((phi - m_phimin) * static_cast<double>(m_invphi));
  float const gphi = 1.0F - fphi;
  // interpolate field values in z, r, phi
  std::array<float, 3> Bzrphi {};
  for (std::size_t i = 0; i < 3; i++) {  // z, r, phi components
    const auto& field = m_field[i];
    Bzrphi[i] = m_scale
        * (gz
               * (gr * (gphi * field[0] + fphi * field[1])
                  + fr * (gphi * field[2] + fphi * field[3]))
           + fz
               * (gr * (gphi * field[4] + fphi * field[5])
                  + fr * (gphi * field[6] + fphi * field[7])));
  }
  // convert (Bz,Br,Bphi) to (Bx,By,Bz)
  const bool positive_r = r > 0.0;
  const float invr = positive_r ? static_cast<float>(1.0 / r) : 0.0F;
  const float c = positive_r ? static_cast<float>(x * static_cast<double>(invr))
                             : static_cast<float>(cos(m_phimin));
  const float es = positive_r
      ? static_cast<float>(y * static_cast<double>(invr))
      : static_cast<float>(sin(m_phimin));
  B[0] = static_cast<double>(Bzrphi[1] * c - Bzrphi[2] * es);
  B[1] = static_cast<double>(Bzrphi[1] * es + Bzrphi[2] * c);
  B[2] = static_cast<double>(Bzrphi[0]);

  // compute field derivatives if requested
  if (deriv != nullptr) {
    float const sz = m_scale * m_invz;
    float const esr = m_scale * m_invr;
    float const sphi = m_scale * m_invphi;
    std::array<float, 3> dBdz {};
    std::array<float, 3> dBdr {};
    std::array<float, 3> dBdphi {};
    for (std::size_t j = 0; j < 3; j++) {  // Bz, Br, Bphi components
      const auto& field = m_field[j];
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
    float const cc = c * c;
    float const cs = c * es;
    float const ss = es * es;
    float const ccinvr = cc * invr;
    float const csinvr = cs * invr;
    float const ssinvr = ss * invr;
    float const sinvr = es * invr;
    float const cinvr = c * invr;
    // the float terms are summed in float, the B[] terms are added in double
    deriv[0] = static_cast<double>(cc * dBdr[1] - cs * dBdr[2]
                                   - csinvr * dBdphi[1] + ssinvr * dBdphi[2])
        + static_cast<double>(sinvr) * B[1];
    deriv[1] = static_cast<double>(cs * dBdr[1] - ss * dBdr[2]
                                   + ccinvr * dBdphi[1] - csinvr * dBdphi[2])
        - static_cast<double>(cinvr) * B[1];
    deriv[2] = static_cast<double>(c * dBdz[1] - es * dBdz[2]);
    deriv[3] = static_cast<double>(cs * dBdr[1] + cc * dBdr[2]
                                   - ssinvr * dBdphi[1] - csinvr * dBdphi[2])
        - static_cast<double>(sinvr) * B[0];
    deriv[4] = static_cast<double>(ss * dBdr[1] + cs * dBdr[2]
                                   + csinvr * dBdphi[1] + ccinvr * dBdphi[2])
        + static_cast<double>(cinvr) * B[0];
    deriv[5] = static_cast<double>(es * dBdz[1] + c * dBdz[2]);
    deriv[6] = static_cast<double>(c * dBdr[0] - sinvr * dBdphi[0]);
    deriv[7] = static_cast<double>(es * dBdr[0] + cinvr * dBdphi[0]);
    deriv[8] = static_cast<double>(dBdz[0]);
  }
}

#endif
