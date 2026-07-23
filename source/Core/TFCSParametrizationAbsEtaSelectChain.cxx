// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#include <cmath>

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
  return val_to_bin(std::abs(extrapol->IDCaloBoundary_eta()));
}

auto TFCSParametrizationAbsEtaSelectChain::get_bin_text(int bin) const
    -> const std::string
{
  if (bin == -1 || bin >= (int)get_number_of_bins()) {
    return std::string(Form("bin=%d not in [%2.2f<=|eta|<%2.2f)",
                            bin,
                            m_bin_low_edge[0],
                            m_bin_low_edge[get_number_of_bins()]));
  }
  return std::string(Form("bin=%d, %2.2f<=|eta|<%2.2f",
                          bin,
                          m_bin_low_edge[bin],
                          m_bin_low_edge[bin + 1]));
}
