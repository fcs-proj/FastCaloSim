/*
  Copyright (C) 2002-2018 CERN for the benefit of the ATLAS collaboration
*/

#include <iostream>

#include "FastCaloSim/TFCSParametrizationEkinSelectChain.h"

#include "CLHEP/Random/RandFlat.h"
#include "FastCaloSim/TFCSExtrapolationState.h"
#include "FastCaloSim/TFCSInvisibleParametrization.h"
#include "FastCaloSim/TFCSSimulationState.h"
#include "FastCaloSim/TFCSTruthState.h"

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

void TFCSParametrizationEkinSelectChain::unit_test(
    TFCSSimulationState* simulstate,
    TFCSTruthState* truth,
    const TFCSExtrapolationState* extrapol)
{
  ISF_FCS::MLogging logger;
  if (!simulstate)
    simulstate = new TFCSSimulationState();
  if (!truth)
    truth = new TFCSTruthState();
  if (!extrapol)
    extrapol = new TFCSExtrapolationState();

  TFCSParametrizationEkinSelectChain chain("chain", "chain");
  chain.setLevel(MSG::DEBUG);

  TFCSParametrization* param;
  param = new TFCSInvisibleParametrization("A begin all", "A begin all");
  param->setLevel(MSG::DEBUG);
  param->set_Ekin_nominal(2);
  param->set_Ekin_min(2);
  param->set_Ekin_max(5);
  chain.push_before_first_bin(param);
  param = new TFCSInvisibleParametrization("A end all", "A end all");
  param->setLevel(MSG::DEBUG);
  param->set_Ekin_nominal(2);
  param->set_Ekin_min(2);
  param->set_Ekin_max(5);
  chain.push_back(param);

  const int n_params = 5;
  for (int i = 2; i < n_params; ++i) {
    param = new TFCSInvisibleParametrization(Form("A%d", i), Form("A %d", i));
    param->setLevel(MSG::DEBUG);
    param->set_Ekin_nominal(TMath::Power(2.0, i));
    param->set_Ekin_min(TMath::Power(2.0, i - 0.5));
    param->set_Ekin_max(TMath::Power(2.0, i + 0.5));
    chain.push_back_in_bin(param);
  }
  for (int i = n_params; i >= 1; --i) {
    param = new TFCSInvisibleParametrization(Form("B%d", i), Form("B %d", i));
    param->setLevel(MSG::DEBUG);
    param->set_Ekin_nominal(TMath::Power(2.0, i));
    param->set_Ekin_min(TMath::Power(2.0, i - 0.5));
    param->set_Ekin_max(TMath::Power(2.0, i + 0.5));
    chain.push_back_in_bin(param);
  }

  ATH_MSG_NOCLASS(logger, "====         Chain setup       ====");
  chain.Print();

  param = new TFCSInvisibleParametrization("B end all", "B end all");
  param->setLevel(MSG::DEBUG);
  chain.push_back(param);
  param = new TFCSInvisibleParametrization("B begin all", "B begin all");
  param->setLevel(MSG::DEBUG);
  chain.push_before_first_bin(param);

  ATH_MSG_NOCLASS(logger, "====         Chain setup       ====");
  chain.Print();
  ATH_MSG_NOCLASS(logger, "==== Simulate with E=0.3      ====");
  truth->SetPtEtaPhiM(0.3, 0, 0, 0);
  chain.simulate(*simulstate, truth, extrapol);
  for (double E = 1; E < 10.1; E += 1) {
    ATH_MSG_NOCLASS(logger, "==== Simulate with E=" << E << "      ====");
    truth->SetPtEtaPhiM(E, 0, 0, 0);
    chain.simulate(*simulstate, truth, extrapol);
  }
  ATH_MSG_NOCLASS(logger, "==== Simulate with E=100      ====");
  truth->SetPtEtaPhiM(100, 0, 0, 0);
  chain.simulate(*simulstate, truth, extrapol);
  ATH_MSG_NOCLASS(logger, "===================================" << std::endl);
  ATH_MSG_NOCLASS(logger, "====== now with random bin ========" << std::endl);
  chain.set_DoRandomInterpolation();
  for (double E = 15; E < 35.1; E += 4) {
    ATH_MSG_NOCLASS(logger, "==== Simulate with E=" << E << "      ====");
    truth->SetPtEtaPhiM(E, 0, 0, 0);
    for (int i = 0; i < 10; ++i)
      chain.simulate(*simulstate, truth, extrapol);
  }
  ATH_MSG_NOCLASS(logger, "===================================" << std::endl);
}
