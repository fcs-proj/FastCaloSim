/*
  Copyright (C) 2002-2018 CERN for the benefit of the ATLAS collaboration
*/

#include "FastCaloSim/Core/TFCSParametrizationEbinChain.h"

#include "FastCaloSim/Core/TFCSExtrapolationState.h"
#include "FastCaloSim/Core/TFCSSimulationState.h"
#include "FastCaloSim/Core/TFCSTruthState.h"

//=============================================
//======= TFCSParametrizationEbinChain =========
//=============================================

const std::string TFCSParametrizationEbinChain::get_variable_text(
    TFCSSimulationState& simulstate,
    const TFCSTruthState*,
    const TFCSExtrapolationState*) const
{
  return std::string(Form("Ebin=%d", simulstate.Ebin()));
}
