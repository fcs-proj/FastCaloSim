// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#include "FastCaloSim/Core/TFCSInvisibleParametrization.h"

//=============================================
//======= TFCSInvisibleParametrization =========
//=============================================

auto TFCSInvisibleParametrization::simulate(
    TFCSSimulationState& /*simulstate*/,
    const TFCSTruthState* /*truth*/,
    const TFCSExtrapolationState* /*extrapol*/) const -> FCSReturnCode
{
  FCS_MSG_VERBOSE(
      "now in TFCSInvisibleParametrization::simulate(). Don't do "
      "anything for invisible");
  return FCSSuccess;
}
