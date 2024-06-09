/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

//
// BFieldCond.h
//
// Current conductors used in BFieldZone.
// It may be finite : p1 to p2.
// It may be infinite : passing p1 and direction parallel to p2.
//
// Masahiro Morii, Harvard University
//
#ifndef BFIELDCOND_H
#define BFIELDCOND_H

#include <cmath>

#include <Eigen/Dense>

class BFieldCond
{
public:
  // constructor
  BFieldCond(bool finite, const double* p1, const double* p2, double curr)
      : m_finite(finite)
      , m_p1(Eigen::Map<const Eigen::Vector3d>(p1))
      , m_p2(Eigen::Map<const Eigen::Vector3d>(p2))
      , m_u(finite ? (m_p2 - m_p1).normalized() : m_p2)
      , m_curr(curr)
      , m_nomCurr(curr)
  {
  }

  // compute magnetic field, plus derivatives if requested, and add
  inline void addBiotSavart(const double* xyz,
                            double* B,
                            double* deriv = 0) const;
  // scale the current wrt the nominal current
  void scaleCurrent(double factor) { m_curr = factor * m_nomCurr; }
  // accessors
  bool finite() const { return m_finite; }
  double p1(int i) const { return m_p1[i]; }
  double p2(int i) const { return m_p2[i]; }
  double curr() const { return m_curr; }

private:
  bool m_finite;  // true if the conductor is finite in length

  /*
   * m_p1: One end of a finite conductor, or one point on an infinite conductor
   * m_p2: The other end of a finite conductor, or the direction vector of an
   * infinite conductor m_u : Direction vector (= m_p2 if infinite)
   */
  const Eigen::Vector3d m_p1, m_p2, m_u;
  double m_curr;  // current (in A) flowing through the conductor
  double m_nomCurr;  // nominal current (in A) read from the map file
};

//
// Compute the Biot-Savart field due to this conductor
// If deriv[9] is passed, also computes the field derivatives
// The results are _added_ to B[] and deriv[].
//
void BFieldCond::addBiotSavart(const double* xyz,
                               double* B_in,
                               double* deriv) const
{
  static const double mu04pi = 1.0e-7;  // mu_0/4pi
  static const double minvsq = 10. * 10.;  // (1 cm)^2

  const Eigen::Map<const Eigen::Vector3d> pos(xyz);
  Eigen::Map<Eigen::Vector3d> B(B_in);

  const Eigen::Vector3d r1 = pos - m_p1;
  const Eigen::Vector3d r2 = pos - m_p2;
  const Eigen::Vector3d v = m_u.cross(r1);
  const double vsq = std::max(v.squaredNorm(), minvsq);

  const double r1u = r1.dot(m_u), r2u = r2.dot(m_u);
  const double r1n = r1.norm(), r2n = r2.norm();

  const double f0 = mu04pi * m_curr / vsq;
  const double f1 = f0 * (m_finite ? r1u / r1n - r2u / r2n : 2.0);
  const double f2 = 2.0 * f1 / vsq;

  B += f1 * v;

  if (deriv) {
    Eigen::Map<Eigen::Matrix<double, 3, 3, Eigen::RowMajor>> D(deriv);

    D(0, 1) -= f1 * m_u(2);
    D(0, 2) += f1 * m_u(1);
    D(1, 0) += f1 * m_u(2);
    D(1, 2) -= f1 * m_u(0);
    D(2, 0) -= f1 * m_u(1);
    D(2, 1) += f1 * m_u(0);

    if (vsq > minvsq) {
      Eigen::Vector3d w = f2 * (m_u * r1u - r1);

      if (m_finite) {
        // Finite conductor segment
        w += f0
            * ((m_u - r1 * r1u / (r1n * r1n)) / r1n
               - (m_u - r2 * r2u / (r2n * r2n)) / r2n);
      }

      D += v * w.transpose();
    }
  }
}

#endif
