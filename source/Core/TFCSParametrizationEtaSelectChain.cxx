// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#include <format>

#include "FastCaloSim/Core/TFCSParametrizationEtaSelectChain.h"

#include "FastCaloSim/Core/TFCSExtrapolationState.h"
#include "FastCaloSim/Core/TFCSSimulationState.h"
#include "FastCaloSim/Core/TFCSTruthState.h"

//=============================================
//======= TFCSParametrizationEtaSelectChain =========
//=============================================

void TFCSParametrizationEtaSelectChain::recalc()
{
  clear();
  if (size() == 0)
    return;

  recalc_pdgid_intersect();
  recalc_Ekin_intersect();
  recalc_eta_union();

  chain().shrink_to_fit();
}

void TFCSParametrizationEtaSelectChain::push_back_in_bin(
    TFCSParametrizationBase* param)
{
  push_back_in_bin(param, param->eta_min(), param->eta_max());
}

auto TFCSParametrizationEtaSelectChain::get_bin(
    TFCSSimulationState&,
    const TFCSTruthState*,
    const TFCSExtrapolationState* extrapol) const -> int
{
  // The bin edges are single precision, so is the binning variable
  return val_to_bin(static_cast<float>(extrapol->IDCaloBoundary_eta()));
}

auto TFCSParametrizationEtaSelectChain::get_variable_text(
    TFCSSimulationState&,
    const TFCSTruthState*,
    const TFCSExtrapolationState* extrapol) const -> const std::string
{
  return std::format("eta={:2.2f}", extrapol->IDCaloBoundary_eta());
}

auto TFCSParametrizationEtaSelectChain::get_bin_text(int bin) const
    -> const std::string
{
  if (bin == -1 || bin >= static_cast<int>(get_number_of_bins())) {
    return std::format("bin={} not in [{:2.2f}<=eta<{:2.2f})",
                       bin,
                       m_bin_low_edge[0],
                       m_bin_low_edge[get_number_of_bins()]);
  }
  return std::format("bin={}, {:2.2f}<=eta<{:2.2f}",
                     bin,
                     m_bin_low_edge[bin],
                     m_bin_low_edge[bin + 1]);
}
