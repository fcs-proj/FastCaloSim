/*
  Copyright (C) 2002-2019 CERN for the benefit of the ATLAS collaboration
*/

#include "FastCaloSim/TFCSParametrizationPlaceholder.h"

//=============================================
//======= TFCSParametrizationPlaceholder =========
//=============================================

FCSReturnCode TFCSParametrizationPlaceholder::simulate(
    TFCSSimulationState& /*simulstate*/,
    const TFCSTruthState* /*truth*/,
    const TFCSExtrapolationState* /*extrapol*/) const
{
  ATH_MSG_ERROR(
      "TFCSParametrizationPlaceholder::simulate(). This is a "
      "placeholder and should never get called. Likely a problem in "
      "the reading of the parametrization file occurred and this "
      "class was not replaced with the real parametrization");
  return FCSFatal;
}
