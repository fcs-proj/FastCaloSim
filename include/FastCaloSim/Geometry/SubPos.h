// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#pragma once

#include <cstdint>

namespace FastCaloSim
{
/// A position within a calorimeter layer along the particle direction. The
/// values are the indices of the sub-positions, e.g. in the extrapolation.
enum class SubPos : std::uint8_t
{
  MID = 0,  // middle
  ENT = 1,  // entrance
  EXT = 2  // exit
};
}  // namespace FastCaloSim
