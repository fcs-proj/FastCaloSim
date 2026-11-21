// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#include <format>

#include "FastCaloSim/Core/TFCSParametrizationEbinChain.h"

#include "FastCaloSim/Core/TFCSExtrapolationState.h"
#include "FastCaloSim/Core/TFCSSimulationState.h"
#include "FastCaloSim/Core/TFCSTruthState.h"

//=============================================
//======= TFCSParametrizationEbinChain =========
//=============================================

auto TFCSParametrizationEbinChain::get_variable_text(
    TFCSSimulationState& simulstate,
    const TFCSTruthState* /*unused*/,
    const TFCSExtrapolationState* /*unused*/) const -> const std::string
{
  return std::format("Ebin={}", bin_label(simulstate.Ebin()));
}
