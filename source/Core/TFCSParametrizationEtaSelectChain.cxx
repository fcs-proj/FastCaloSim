/*
  Copyright (C) 2002-2018 CERN for the benefit of the ATLAS collaboration
*/

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

int TFCSParametrizationEtaSelectChain::get_bin(
    TFCSSimulationState&,
    const TFCSTruthState*,
    const TFCSExtrapolationState* extrapol) const
{
  return val_to_bin(extrapol->IDCaloBoundary_eta());
}

const std::string TFCSParametrizationEtaSelectChain::get_variable_text(
    TFCSSimulationState&,
    const TFCSTruthState*,
    const TFCSExtrapolationState* extrapol) const
{
  return std::string(Form("eta=%2.2f", extrapol->IDCaloBoundary_eta()));
}

const std::string TFCSParametrizationEtaSelectChain::get_bin_text(int bin) const
{
  if (bin == -1 || bin >= (int)get_number_of_bins()) {
    return std::string(Form("bin=%d not in [%2.2f<=eta<%2.2f)",
                            bin,
                            m_bin_low_edge[0],
                            m_bin_low_edge[get_number_of_bins()]));
  }
  return std::string(Form("bin=%d, %2.2f<=eta<%2.2f",
                          bin,
                          m_bin_low_edge[bin],
                          m_bin_low_edge[bin + 1]));
}
