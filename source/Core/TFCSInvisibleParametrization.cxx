// Copyright (c) 2024 CERN for the benefit of the FastCaloSim project

#include "FastCaloSim/Core/TFCSInvisibleParametrization.h"

//=============================================
//======= TFCSInvisibleParametrization =========
//=============================================

FCSReturnCode TFCSInvisibleParametrization::simulate(
    TFCSSimulationState& /*simulstate*/,
    const TFCSTruthState* /*truth*/,
    const TFCSExtrapolationState* /*extrapol*/) const
{
  FCS_MSG_VERBOSE(
      "now in TFCSInvisibleParametrization::simulate(). Don't do "
      "anything for invisible");
  return FCSSuccess;
}
