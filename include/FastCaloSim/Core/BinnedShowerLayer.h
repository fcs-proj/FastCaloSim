// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#pragma once

#include <vector>

namespace FastCaloSim
{
/// One calorimeter layer of a binned shower: the indices of the bins with an
/// energy deposit and the fraction of the shower energy in each of them
struct BinnedShowerLayer
{
  std::vector<unsigned int> bin_index_vector;
  std::vector<float> E_vector;
};
}  // namespace FastCaloSim
