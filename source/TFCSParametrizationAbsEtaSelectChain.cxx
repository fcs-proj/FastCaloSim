/*
  Copyright (C) 2002-2018 CERN for the benefit of the ATLAS collaboration
*/

#include <iostream>

#include "FastCaloSim/TFCSParametrizationAbsEtaSelectChain.h"

#include "FastCaloSim/TFCSExtrapolationState.h"
#include "FastCaloSim/TFCSSimulationState.h"
#include "FastCaloSim/TFCSTruthState.h"

//=============================================
//======= TFCSParametrizationAbsEtaSelectChain =========
//=============================================

int TFCSParametrizationAbsEtaSelectChain::get_bin(
    TFCSSimulationState&,
    const TFCSTruthState*,
    const TFCSExtrapolationState* extrapol) const
{
  return val_to_bin(TMath::Abs(extrapol->IDCaloBoundary_eta()));
}

const std::string TFCSParametrizationAbsEtaSelectChain::get_bin_text(
    int bin) const
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

void TFCSParametrizationAbsEtaSelectChain::unit_test(
    TFCSSimulationState* simulstate,
    TFCSTruthState* truth,
    TFCSExtrapolationState* extrapol)
{
  ISF_FCS::MLogging logger;
  if (!simulstate)
    simulstate = new TFCSSimulationState();
  if (!truth)
    truth = new TFCSTruthState();
  if (!extrapol)
    extrapol = new TFCSExtrapolationState();

  TFCSParametrizationAbsEtaSelectChain chain("chain", "chain");
  chain.setLevel(MSG::DEBUG);

  TFCSParametrization* param;
  param = new TFCSParametrization("A begin all", "A begin all");
  param->setLevel(MSG::DEBUG);
  param->set_eta_nominal(2);
  param->set_eta_min(2);
  param->set_eta_max(3);
  chain.push_before_first_bin(param);
  param = new TFCSParametrization("A end all", "A end all");
  param->setLevel(MSG::DEBUG);
  param->set_eta_nominal(2);
  param->set_eta_min(2);
  param->set_eta_max(3);
  chain.push_back(param);

  const int n_params = 5;
  for (int i = 2; i < n_params; ++i) {
    param = new TFCSParametrization(Form("A%d", i), Form("A %d", i));
    param->setLevel(MSG::DEBUG);
    param->set_eta_nominal(i * i + 0.1);
    param->set_eta_min(i * i);
    param->set_eta_max((i + 1) * (i + 1));
    chain.push_back_in_bin(param);
  }
  for (int i = n_params; i >= 1; --i) {
    param = new TFCSParametrization(Form("B%d", i), Form("B %d", i));
    param->setLevel(MSG::DEBUG);
    param->set_eta_nominal(i * i + 0.1);
    param->set_eta_min(i * i);
    param->set_eta_max((i + 1) * (i + 1));
    chain.push_back_in_bin(param);
  }

  ATH_MSG_NOCLASS(logger, "====         Chain setup       ====");
  chain.Print();

  param = new TFCSParametrization("B end all", "B end all");
  param->setLevel(MSG::DEBUG);
  chain.push_back(param);
  param = new TFCSParametrization("B begin all", "B begin all");
  param->setLevel(MSG::DEBUG);
  chain.push_before_first_bin(param);

  ATH_MSG_NOCLASS(logger, "====         Chain setup       ====");
  chain.Print();
  ATH_MSG_NOCLASS(logger, "==== Simulate with eta=0.1      ====");
  extrapol->set_IDCaloBoundary_eta(0.1);
  chain.simulate(*simulstate, truth, extrapol);
  ATH_MSG_NOCLASS(logger, "==== Simulate with eta=-1.1      ====");
  extrapol->set_IDCaloBoundary_eta(-1.1);
  chain.simulate(*simulstate, truth, extrapol);
  ATH_MSG_NOCLASS(logger, "==== Simulate with eta=2.1      ====");
  extrapol->set_IDCaloBoundary_eta(2.1);
  chain.simulate(*simulstate, truth, extrapol);
  ATH_MSG_NOCLASS(logger, "==== Simulate with eta=4.1      ====");
  extrapol->set_IDCaloBoundary_eta(4.1);
  chain.simulate(*simulstate, truth, extrapol);
  ATH_MSG_NOCLASS(logger, "==== Simulate with eta=-4.1      ====");
  extrapol->set_IDCaloBoundary_eta(-4.1);
  chain.simulate(*simulstate, truth, extrapol);
  ATH_MSG_NOCLASS(logger, "==== Simulate with eta=100      ====");
  extrapol->set_IDCaloBoundary_eta(100);
  chain.simulate(*simulstate, truth, extrapol);
  ATH_MSG_NOCLASS(logger, "===================================" << std::endl);
}
