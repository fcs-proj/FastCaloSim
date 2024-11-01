// Copyright (c) 2024 CERN for the benefit of the FastCaloSim project

#include "FastCaloSim/Core/TFCSInitWithEkin.h"

#include "FastCaloSim/Core/TFCSSimulationState.h"
#include "FastCaloSim/Core/TFCSTruthState.h"

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
  FCS_MSG_DEBUG("set E to Ekin=" << truth->Ekin());
  simulstate.set_E(truth->Ekin());
  return FCSSuccess;
}
