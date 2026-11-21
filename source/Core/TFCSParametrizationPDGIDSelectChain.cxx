// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#include <iostream>

#include "FastCaloSim/Core/TFCSParametrizationPDGIDSelectChain.h"

#include "FastCaloSim/Core/TFCSExtrapolationState.h"
#include "FastCaloSim/Core/TFCSSimulationState.h"
#include "FastCaloSim/Core/TFCSTruthState.h"

//=============================================
//======= TFCSParametrizationPDGIDSelectChain =========
//=============================================

void TFCSParametrizationPDGIDSelectChain::recalc()
{
  clear();
  if (size() == 0)
    return;

  recalc_pdgid_union();
  recalc_Ekin_eta_intersect();

  chain().shrink_to_fit();
}

auto TFCSParametrizationPDGIDSelectChain::simulate(
    TFCSSimulationState& simulstate,
    const TFCSTruthState* truth,
    const TFCSExtrapolationState* extrapol) const -> FCSReturnCode
{
  return simulate_with_restarts(
      "TFCSParametrizationPDGIDSelectChain::simulate()",
      [&]
      {
        FCS_MSG_DEBUG("Running for pdgid=" << truth->pdgid());
        for (const auto& param : chain()) {
          FCS_MSG_DEBUG("Now testing: "
                        << param->GetName()
                        << ((SimulateOnlyOnePDGID() == true)
                                ? ", abort PDGID loop afterwards"
                                : ", continue PDGID loop afterwards"));
          if (!param->is_match_pdgid(truth->pdgid()))
            continue;
          FCS_MSG_DEBUG("pdgid=" << truth->pdgid()
                                 << ", now run: " << param->GetName()
                                 << ((SimulateOnlyOnePDGID() == true)
                                         ? ", abort PDGID loop afterwards"
                                         : ", continue PDGID loop afterwards"));
          const FCSReturnCode status =
              simulate_and_retry(param, simulstate, truth, extrapol);
          if (status != FCSSuccess)
            return status;
          if (SimulateOnlyOnePDGID())
            break;
        }
        return FCSSuccess;
      });
}
