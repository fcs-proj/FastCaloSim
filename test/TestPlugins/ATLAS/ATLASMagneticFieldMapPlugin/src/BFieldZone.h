/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

//
// BFieldZone.h
//
// A "zone" inside the toroid field map
//
// Masahiro Morii, Harvard University
//
#ifndef BFIELDZONE_H
#define BFIELDZONE_H

#include <vector>

#include "BFieldCond.h"
#include "BFieldMesh.h"

class BFieldZone : public BFieldMesh<short>
{
public:
  // constructor
  BFieldZone(int id,
             double zmin,
             double zmax,
             double rmin,
             double rmax,
             double phimin,
             double phimax,
             double scale)
      : BFieldMesh<short>(zmin, zmax, rmin, rmax, phimin, phimax, scale)
      , m_id(id)
  {
    ;
  }
  // add elements to vectors
  void appendCond(const BFieldCond& cond) { m_cond.push_back(cond); }
  // compute Biot-Savart magnetic field and add to B[3]
  inline void addBiotSavart(const double* xyz,
                            double* B,
                            double* deriv = nullptr) const;
  // scale B field by a multiplicative factor
  void scaleField(double factor)
  {
    scaleBscale(factor);
    for (unsigned i = 0; i < ncond(); i++) {
      m_cond[i].scaleCurrent(factor);
    }
  }
  // accessors
  auto id() const -> int { return m_id; }
  auto ncond() const -> unsigned
  {
    return static_cast<unsigned>(m_cond.size());
  }
  auto cond(int i) const -> const BFieldCond&
  {
    return m_cond[static_cast<std::size_t>(i)];
  }
  auto condVector() const -> const std::vector<BFieldCond>* { return &m_cond; }
  auto memSize() const -> int
  {
    return static_cast<int>(
        static_cast<std::size_t>(BFieldMesh<short>::memSize()) + sizeof(int)
        + sizeof(BFieldCond) * m_cond.capacity());
  }
  // adjust the min/max edges to a new value
  void adjustMin(int i, double x)
  {
    m_min[i] = x;
    m_mesh[i].front() = x;
  }
  void adjustMax(int i, double x)
  {
    m_max[i] = x;
    m_mesh[i].back() = x;
  }

private:
  int m_id;  // zone ID number
  std::vector<BFieldCond> m_cond;  // list of current conductors
};

// inline functions

//
// Compute magnetic field due to the conductors
//
void BFieldZone::addBiotSavart(const double* xyz,
                               double* B,
                               double* deriv) const
{
  for (const auto& i : m_cond) {
    i.addBiotSavart(xyz, B, deriv);
  }
}

#endif
