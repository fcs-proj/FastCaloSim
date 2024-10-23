// Copyright (c) 2024 CERN for the benefit of the FastCaloSim project

#include <iostream>

#include "FastCaloSim/Core/TFCSParametrizationEkinSelectChain.h"

#include "CLHEP/Random/RandFlat.h"
#include "FastCaloSim/Core/TFCSExtrapolationState.h"
#include "FastCaloSim/Core/TFCSSimulationState.h"
#include "FastCaloSim/Core/TFCSTruthState.h"

//=============================================
//======= TFCSParametrizationEkinSelectChain =========
//=============================================

void TFCSParametrizationEkinSelectChain::recalc()
{
  clear();
  if (size() == 0)
    return;

  recalc_pdgid_intersect();
  recalc_Ekin_union();
  recalc_eta_intersect();

  chain().shrink_to_fit();
}

void TFCSParametrizationEkinSelectChain::push_back_in_bin(
    TFCSParametrizationBase* param)
{
  push_back_in_bin(param, param->Ekin_min(), param->Ekin_max());
}

int TFCSParametrizationEkinSelectChain::get_bin(
    TFCSSimulationState& simulstate,
    const TFCSTruthState* truth,
    const TFCSExtrapolationState*) const
{
  if (!simulstate.randomEngine()) {
    return -1;
  }

  float Ekin = truth->Ekin();
  int bin = val_to_bin(Ekin);

  if (!DoRandomInterpolation())
    return bin;
  float rnd = CLHEP::RandFlat::shoot(simulstate.randomEngine());
  if (bin < 0)
    return bin;
  if (bin >= (int)get_number_of_bins())
    return bin;

  // if no parametrizations for this bin, return
  if (m_bin_start[bin + 1] == m_bin_start[bin])
    return bin;

  TFCSParametrizationBase* first_in_bin = chain()[m_bin_start[bin]];
  if (!first_in_bin)
    return bin;

  if (Ekin < first_in_bin->Ekin_nominal()) {
    if (bin == 0)
      return bin;
    int prevbin = bin - 1;
    // if no parametrizations for previous bin, return
    if (m_bin_start[prevbin + 1] == m_bin_start[prevbin])
      return bin;

    TFCSParametrizationBase* first_in_prevbin = chain()[m_bin_start[prevbin]];
    if (!first_in_prevbin)
      return bin;

    float logEkin = TMath::Log(Ekin);
    float logEkin_nominal = TMath::Log(first_in_bin->Ekin_nominal());
    float logEkin_previous = TMath::Log(first_in_prevbin->Ekin_nominal());
    float numerator = logEkin - logEkin_previous;
    float denominator = logEkin_nominal - logEkin_previous;
    if (denominator <= 0)
      return bin;

    if (numerator / denominator < rnd)
      bin = prevbin;
    ATH_MSG_DEBUG("logEkin="
                  << logEkin << " logEkin_previous=" << logEkin_previous
                  << " logEkin_nominal=" << logEkin_nominal
                  << " (rnd=" << 1 - rnd
                  << " < p(previous)=" << (1 - numerator / denominator)
                  << ")? => orgbin=" << prevbin + 1 << " selbin=" << bin);
  } else {
    if (bin == (int)get_number_of_bins() - 1)
      return bin;
    int nextbin = bin + 1;
    // if no parametrizations for previous bin, return
    if (m_bin_start[nextbin + 1] == m_bin_start[nextbin])
      return bin;

    TFCSParametrizationBase* first_in_nextbin = chain()[m_bin_start[nextbin]];
    if (!first_in_nextbin)
      return bin;

    float logEkin = TMath::Log(Ekin);
    float logEkin_nominal = TMath::Log(first_in_bin->Ekin_nominal());
    float logEkin_next = TMath::Log(first_in_nextbin->Ekin_nominal());
    float numerator = logEkin - logEkin_nominal;
    float denominator = logEkin_next - logEkin_nominal;
    if (denominator <= 0)
      return bin;

    if (rnd < numerator / denominator)
      bin = nextbin;
    ATH_MSG_DEBUG("logEkin="
                  << logEkin << " logEkin_nominal=" << logEkin_nominal
                  << " logEkin_next=" << logEkin_next << " (rnd=" << rnd
                  << " < p(next)=" << numerator / denominator
                  << ")? => orgbin=" << nextbin - 1 << " selbin=" << bin);
  }

  return bin;
}

const std::string TFCSParametrizationEkinSelectChain::get_variable_text(
    TFCSSimulationState&,
    const TFCSTruthState* truth,
    const TFCSExtrapolationState*) const
{
  return std::string(Form("Ekin=%1.1f", truth->Ekin()));
}

const std::string TFCSParametrizationEkinSelectChain::get_bin_text(
    int bin) const
{
  if (bin == -1 || bin >= (int)get_number_of_bins()) {
    return std::string(Form("bin=%d not in [%1.1f<=Ekin<%1.1f)",
                            bin,
                            m_bin_low_edge[0],
                            m_bin_low_edge[get_number_of_bins()]));
  }
  if (DoRandomInterpolation()) {
    return std::string(Form("bin=%d, %1.1f<=Ekin(+random)<%1.1f",
                            bin,
                            m_bin_low_edge[bin],
                            m_bin_low_edge[bin + 1]));
  }
  return std::string(Form("bin=%d, %1.1f<=Ekin<%1.1f",
                          bin,
                          m_bin_low_edge[bin],
                          m_bin_low_edge[bin + 1]));
}
