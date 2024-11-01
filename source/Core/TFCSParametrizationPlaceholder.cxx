// Copyright (c) 2024 CERN for the benefit of the FastCaloSim project

#include "FastCaloSim/Core/TFCSParametrizationPlaceholder.h"

//=============================================
//======= TFCSParametrizationPlaceholder =========
//=============================================

FCSReturnCode TFCSParametrizationPlaceholder::simulate(
    TFCSSimulationState& /*simulstate*/,
    const TFCSTruthState* /*truth*/,
    const TFCSExtrapolationState* /*extrapol*/) const
{
  FCS_MSG_ERROR(
      "TFCSParametrizationPlaceholder::simulate(). This is a "
      "placeholder and should never get called. Likely a problem in "
      "the reading of the parametrization file occurred and this "
      "class was not replaced with the real parametrization");
  return FCSFatal;
}
