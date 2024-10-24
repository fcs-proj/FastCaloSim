// Copyright (c) 2024 CERN for the benefit of the FastCaloSim project

/*
  Class to store particle charge data. As not all possible particles are defined
  here, experiments must provide their own implementation to ensure full
  particle coverage. A future solution will transition to more generic
  approaches, such as those in the AtlasPID header:
  https://gitlab.cern.ch/atlas/athena/-/blob/main/Generators/TruthUtils/TruthUtils/AtlasPID.h?ref_type=heads
*/

#pragma once

#include <stdexcept>

namespace ParticleData
{

constexpr auto charge(const int pdgID) -> double
{
  if (pdgID == 11 || pdgID == 211 || pdgID == 2212)
    return 1.;
  else if (pdgID == -11 || pdgID == -211 || pdgID == -2212)
    return -1.;
  else if (pdgID == 22 || std::abs(pdgID) == 2112 || pdgID == 111)
    return 0;
  else {
    std::runtime_error("Error: This pdgID is not supported: ");
  }
}

}  // namespace ParticleData
