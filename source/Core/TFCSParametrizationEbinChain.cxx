// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#include "FastCaloSim/Core/TFCSParametrizationEbinChain.h"

#include "FastCaloSim/Core/TFCSExtrapolationState.h"
#include "FastCaloSim/Core/TFCSSimulationState.h"
#include "FastCaloSim/Core/TFCSTruthState.h"

//=============================================
//======= TFCSParametrizationEbinChain =========
//=============================================

auto TFCSParametrizationEbinChain::get_variable_text(
    TFCSSimulationState& simulstate,
    const TFCSTruthState*,
    const TFCSExtrapolationState*) const -> const std::string
{
  return std::string(Form("Ebin=%d", simulstate.Ebin()));
}
