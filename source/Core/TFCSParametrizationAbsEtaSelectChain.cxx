// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#include <cmath>
#include <format>

#include "FastCaloSim/Core/TFCSParametrizationAbsEtaSelectChain.h"

#include "FastCaloSim/Core/TFCSExtrapolationState.h"
#include "FastCaloSim/Core/TFCSSimulationState.h"
#include "FastCaloSim/Core/TFCSTruthState.h"

//=============================================
//======= TFCSParametrizationAbsEtaSelectChain =========
//=============================================

auto TFCSParametrizationAbsEtaSelectChain::get_bin(
    TFCSSimulationState&,
    const TFCSTruthState*,
    const TFCSExtrapolationState* extrapol) const -> int
{
  // The bin edges are single precision, so is the binning variable
  return val_to_bin(
      static_cast<float>(std::abs(extrapol->IDCaloBoundary_eta())));
}

auto TFCSParametrizationAbsEtaSelectChain::get_bin_text(int bin) const
    -> const std::string
{
  if (bin == -1 || bin >= static_cast<int>(get_number_of_bins())) {
    return std::format("bin={} not in [{:2.2f}<=|eta|<{:2.2f})",
                       bin,
                       m_bin_low_edge[0],
                       m_bin_low_edge[get_number_of_bins()]);
  }
  return std::format("bin={}, {:2.2f}<=|eta|<{:2.2f}",
                     bin,
                     m_bin_low_edge[bin],
                     m_bin_low_edge[bin + 1]);
}
