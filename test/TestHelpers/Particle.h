// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#pragma once

#include "G4ThreeVector.hh"

namespace TestHelpers
{

class RPhiEta
{
public:
  RPhiEta(float r, float phi, float eta)
  {
    m_dir.setRhoPhiEta(static_cast<double>(r),
                       static_cast<double>(phi),
                       static_cast<double>(eta));
  }
  inline auto get() const -> const G4ThreeVector& { return m_dir; }

private:
  G4ThreeVector m_dir;
};

class RZPhi
{
public:
  RZPhi(float r, float z, float phi)
  {
    m_dir.setRhoPhiZ(static_cast<double>(r),
                     static_cast<double>(phi),
                     static_cast<double>(z));
  }
  inline auto get() const -> const G4ThreeVector& { return m_dir; }

private:
  G4ThreeVector m_dir;
};

struct Particle
{
  int pid;
  float ekin;
  G4ThreeVector pos;
  G4ThreeVector dir;
};

}  // namespace TestHelpers
