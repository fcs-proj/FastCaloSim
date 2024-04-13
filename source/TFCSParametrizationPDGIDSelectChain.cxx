/*
  Copyright (C) 2002-2020 CERN for the benefit of the ATLAS collaboration
*/

#include "FastCaloSim/TFCSParametrizationPDGIDSelectChain.h"
#include "FastCaloSim/TFCSInvisibleParametrization.h"
#include "FastCaloSim/TFCSSimulationState.h"
#include "FastCaloSim/TFCSTruthState.h"
#include "FastCaloSim/TFCSExtrapolationState.h"
#include <iostream>

//=============================================
//======= TFCSParametrizationPDGIDSelectChain =========
//=============================================

void TFCSParametrizationPDGIDSelectChain::recalc() {
  clear();
  if (size() == 0)
    return;

  recalc_pdgid_union();
  recalc_Ekin_eta_intersect();

  chain().shrink_to_fit();
}

FCSReturnCode TFCSParametrizationPDGIDSelectChain::simulate(
    TFCSSimulationState &simulstate, const TFCSTruthState *truth,
    const TFCSExtrapolationState *extrapol) const {
  Int_t retry = 0;
  Int_t retry_warning = 1;

  FCSReturnCode status = FCSSuccess;
  for (int i = 0; i <= retry; i++) {
    if (i >= retry_warning)
      ATH_MSG_WARNING("TFCSParametrizationPDGIDSelectChain::simulate(): Retry "
                      "simulate call "
                      << i << "/" << retry);

    ATH_MSG_DEBUG("Running for pdgid=" << truth->pdgid());
    for (const auto &param : chain()) {
      ATH_MSG_DEBUG("Now testing: "
                    << param->GetName()
                    << ((SimulateOnlyOnePDGID() == true)
                            ? ", abort PDGID loop afterwards"
                            : ", continue PDGID loop afterwards"));
      if (param->is_match_pdgid(truth->pdgid())) {
        ATH_MSG_DEBUG("pdgid=" << truth->pdgid()
                               << ", now run: " << param->GetName()
                               << ((SimulateOnlyOnePDGID() == true)
                                       ? ", abort PDGID loop afterwards"
                                       : ", continue PDGID loop afterwards"));
        status = simulate_and_retry(param, simulstate, truth, extrapol);
        if (status >= FCSRetry) {
          retry = status - FCSRetry;
          retry_warning = retry >> 1;
          if (retry_warning < 1)
            retry_warning = 1;
          break;
        }
        if (status == FCSFatal)
          return FCSFatal;

        if (SimulateOnlyOnePDGID())
          break;
      }
    }

    if (status == FCSSuccess)
      break;
  }

  if (status != FCSSuccess) {
    ATH_MSG_FATAL(
        "TFCSParametrizationChain::simulate(): Simulate call failed after "
        << retry << " retries");
    return FCSFatal;
  }

  return FCSSuccess;
}

void TFCSParametrizationPDGIDSelectChain::unit_test(
    TFCSSimulationState *simulstate, TFCSTruthState *truth,
    TFCSExtrapolationState *extrapol) {
  ISF_FCS::MLogging logger;
  if (!simulstate)
    simulstate = new TFCSSimulationState();
  if (!truth)
    truth = new TFCSTruthState();
  if (!extrapol)
    extrapol = new TFCSExtrapolationState();

  TFCSParametrizationPDGIDSelectChain chain("chain", "chain");
  chain.setLevel(MSG::DEBUG);

  TFCSParametrization *param;
  param = new TFCSInvisibleParametrization("A begin all", "A begin all");
  param->setLevel(MSG::VERBOSE);
  param->set_pdgid(0);
  chain.push_back(param);
  param = new TFCSInvisibleParametrization("A end all", "A end all");
  param->setLevel(MSG::VERBOSE);
  param->set_pdgid(0);
  chain.push_back(param);

  for (int i = 0; i < 3; ++i) {
    param = new TFCSInvisibleParametrization(Form("A%d", i), Form("A %d", i));
    param->setLevel(MSG::VERBOSE);
    param->set_pdgid(i);
    chain.push_back(param);
  }

  for (int i = 3; i > 0; --i) {
    param = new TFCSParametrization(Form("B%d", i), Form("B %d", i));
    param->setLevel(MSG::VERBOSE);
    param->set_pdgid(i);
    chain.push_back(param);
  }
  param = new TFCSInvisibleParametrization("B end all", "B end all");
  param->setLevel(MSG::VERBOSE);
  param->set_match_all_pdgid();
  chain.push_back(param);
  param = new TFCSInvisibleParametrization("B begin all", "B begin all");
  param->setLevel(MSG::VERBOSE);
  param->set_pdgid(1);
  chain.push_back(param);

  ATH_MSG_NOCLASS(logger, "====         Chain setup       ====");
  chain.Print();
  ATH_MSG_NOCLASS(logger, "==== Simulate with pdgid=0      ====");
  truth->set_pdgid(0);
  chain.simulate(*simulstate, truth, extrapol);
  ATH_MSG_NOCLASS(logger, "==== Simulate with pdgid=1      ====");
  truth->set_pdgid(1);
  chain.simulate(*simulstate, truth, extrapol);
  ATH_MSG_NOCLASS(logger, "==== Simulate with pdgid=2      ====");
  truth->set_pdgid(2);
  chain.set_RetryChainFromStart();
  chain.simulate(*simulstate, truth, extrapol);
  chain.reset_RetryChainFromStart();
  ATH_MSG_NOCLASS(logger, "=====================================" << std::endl);

  ATH_MSG_NOCLASS(logger, "=====================================");
  ATH_MSG_NOCLASS(logger, "= Now only one simul for each PDGID =");
  ATH_MSG_NOCLASS(logger, "=====================================");
  chain.set_SimulateOnlyOnePDGID();
  ATH_MSG_NOCLASS(logger, "==== Simulate with pdgid=0      ====");
  truth->set_pdgid(0);
  chain.simulate(*simulstate, truth, extrapol);
  ATH_MSG_NOCLASS(logger, "==== Simulate with pdgid=1      ====");
  truth->set_pdgid(1);
  chain.simulate(*simulstate, truth, extrapol);
  ATH_MSG_NOCLASS(logger, "==== Simulate with pdgid=2      ====");
  truth->set_pdgid(2);
  chain.set_RetryChainFromStart();
  chain.simulate(*simulstate, truth, extrapol);
  chain.reset_RetryChainFromStart();
}
