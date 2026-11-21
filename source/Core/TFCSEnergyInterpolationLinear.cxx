// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

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

{
}

auto TFCSEnergyInterpolationLinear::simulate(
    TFCSSimulationState& simulstate,
    const TFCSTruthState* truth,
    const TFCSExtrapolationState* /*extrapol*/) const -> FCSReturnCode
{
  const double Emean = static_cast<double>(m_slope) * truth->Ekin()
      + static_cast<double>(m_offset);

  FCS_MSG_DEBUG("set E=" << Emean << " for true Ekin=" << truth->Ekin());
  simulstate.set_E(Emean);

  return FCSSuccess;
}

void TFCSEnergyInterpolationLinear::Print(Option_t* option) const
{
  const auto [shortprint, longprint, optprint] = print_options(option);
  TFCSParametrization::Print(option);

  if (longprint) {
    FCS_MSG_INFO(optprint << "  Emean=" << m_slope << "*Ekin(true) + "
                          << m_offset);
  }
}
