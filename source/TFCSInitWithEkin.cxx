/*
  Copyright (C) 2002-2019 CERN for the benefit of the ATLAS collaboration
*/

#include "FastCaloSim/TFCSInitWithEkin.h"

#include "FastCaloSim/TFCSSimulationState.h"
#include "FastCaloSim/TFCSTruthState.h"

//=============================================
//======= TFCSInitWithEkin =========
//=============================================

TFCSInitWithEkin::TFCSInitWithEkin(const char* name, const char* title)
    : TFCSParametrization(name, title)
{
  set_match_all_pdgid();
}

FCSReturnCode TFCSInitWithEkin::simulate(TFCSSimulationState& simulstate,
                                         const TFCSTruthState* truth,
                                         const TFCSExtrapolationState*) const
{
  ATH_MSG_DEBUG("set E to Ekin=" << truth->Ekin());
  simulstate.set_E(truth->Ekin());
  return FCSSuccess;
}
