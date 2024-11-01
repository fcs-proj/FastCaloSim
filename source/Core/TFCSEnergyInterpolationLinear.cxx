// Copyright (c) 2024 CERN for the benefit of the FastCaloSim project

#include <iostream>

#include "FastCaloSim/Core/TFCSEnergyInterpolationLinear.h"

#include "FastCaloSim/Core/TFCSExtrapolationState.h"
#include "FastCaloSim/Core/TFCSSimulationState.h"
#include "FastCaloSim/Core/TFCSTruthState.h"

//=============================================
//======= TFCSEnergyInterpolationLinear =========
//=============================================

TFCSEnergyInterpolationLinear::TFCSEnergyInterpolationLinear(const char* name,
                                                             const char* title)
    : TFCSParametrization(name, title)
    , m_slope(1)
    , m_offset(0)
{
}

FCSReturnCode TFCSEnergyInterpolationLinear::simulate(
    TFCSSimulationState& simulstate,
    const TFCSTruthState* truth,
    const TFCSExtrapolationState*) const
{
  const float Emean = m_slope * truth->Ekin() + m_offset;

  FCS_MSG_DEBUG("set E=" << Emean << " for true Ekin=" << truth->Ekin());
  simulstate.set_E(Emean);

  return FCSSuccess;
}

void TFCSEnergyInterpolationLinear::Print(Option_t* option) const
{
  TString opt(option);
  bool shortprint = opt.Index("short") >= 0;
  bool longprint =
      msgLvl(FCS_MSG::DEBUG) || (msgLvl(FCS_MSG::INFO) && !shortprint);
  TString optprint = opt;
  optprint.ReplaceAll("short", "");
  TFCSParametrization::Print(option);

  if (longprint)
    FCS_MSG_INFO(optprint << "  Emean=" << m_slope << "*Ekin(true) + "
                          << m_offset);
}
