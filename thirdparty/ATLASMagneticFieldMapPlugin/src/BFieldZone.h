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
                            double* deriv = 0) const;
  // scale B field by a multiplicative factor
  void scaleField(double factor)
  {
    scaleBscale(factor);
    for (unsigned i = 0; i < ncond(); i++) {
      m_cond[i].scaleCurrent(factor);
    }
  }
  // accessors
  int id() const { return m_id; }
  unsigned ncond() const { return m_cond.size(); }
  const BFieldCond& cond(int i) const { return m_cond[i]; }
  const std::vector<BFieldCond>* condVector() const { return &m_cond; }
  int memSize() const
  {
    return BFieldMesh<short>::memSize() + sizeof(int)
        + sizeof(BFieldCond) * m_cond.capacity();
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
  for (unsigned i = 0; i < m_cond.size(); i++) {
    m_cond[i].addBiotSavart(xyz, B, deriv);
  }
}

#endif
