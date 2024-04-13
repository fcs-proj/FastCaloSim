/*
  Copyright (C) 2002-2019 CERN for the benefit of the ATLAS collaboration
*/

#include "FastCaloSim/TFCSInvisibleParametrization.h"

//=============================================
//======= TFCSInvisibleParametrization =========
//=============================================

FCSReturnCode TFCSInvisibleParametrization::simulate(
    TFCSSimulationState & /*simulstate*/, const TFCSTruthState * /*truth*/,
    const TFCSExtrapolationState * /*extrapol*/) const {
  ATH_MSG_VERBOSE("now in TFCSInvisibleParametrization::simulate(). Don't do "
                  "anything for invisible");
  return FCSSuccess;
}
