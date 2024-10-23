// Copyright (c) 2024 CERN for the benefit of the FastCaloSim project

#pragma once

#include "G4ThreeVector.hh"

namespace TestHelpers
{

class RPhiEta
{
public:
  RPhiEta(float r, float phi, float eta) { m_dir.setRhoPhiEta(r, phi, eta); }
  inline auto get() const -> const G4ThreeVector& { return m_dir; }

private:
  G4ThreeVector m_dir;
};

class RZPhi
{
public:
  RZPhi(float r, float z, float phi) { m_dir.setRhoPhiZ(r, phi, z); }
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
