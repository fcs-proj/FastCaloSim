// Copyright (c) 2024 CERN for the benefit of the FastCaloSim project

#include <iostream>

#include "FastCaloSim/Core/TFCSParametrizationPDGIDSelectChain.h"

#include "FastCaloSim/Core/TFCSExtrapolationState.h"
#include "FastCaloSim/Core/TFCSSimulationState.h"
#include "FastCaloSim/Core/TFCSTruthState.h"

//=============================================
//======= TFCSParametrizationPDGIDSelectChain =========
//=============================================
using namespace FastCaloSim::Core;

void TFCSParametrizationPDGIDSelectChain::recalc()
{
  clear();
  if (size() == 0)
    return;

  recalc_pdgid_union();
  recalc_Ekin_eta_intersect();

  chain().shrink_to_fit();
}

FCSReturnCode TFCSParametrizationPDGIDSelectChain::simulate(
    TFCSSimulationState& simulstate,
    const TFCSTruthState* truth,
    const TFCSExtrapolationState* extrapol) const
{
  Int_t retry = 0;
  Int_t retry_warning = 1;

  FCSReturnCode status = FCSSuccess;
  for (int i = 0; i <= retry; i++) {
    if (i >= retry_warning)
      MSG_WARNING(
          "TFCSParametrizationPDGIDSelectChain::simulate(): Retry "
          "simulate call "
          << i << "/" << retry);

    MSG_DEBUG("Running for pdgid=" << truth->pdgid());
    for (const auto& param : chain()) {
      MSG_DEBUG("Now testing: " << param->GetName()
                                << ((SimulateOnlyOnePDGID() == true)
                                        ? ", abort PDGID loop afterwards"
                                        : ", continue PDGID loop afterwards"));
      if (param->is_match_pdgid(truth->pdgid())) {
        MSG_DEBUG("pdgid=" << truth->pdgid()
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
    MSG_FATAL(
        "TFCSParametrizationChain::simulate(): Simulate call failed after "
        << retry << " retries");
    return FCSFatal;
  }

  return FCSSuccess;
}
