// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#include <iostream>

#include "FastCaloSim/Core/TFCSEnergyInterpolationHistogram.h"

#include "FastCaloSim/Core/TFCSExtrapolationState.h"
#include "FastCaloSim/Core/TFCSSimulationState.h"
#include "FastCaloSim/Core/TFCSTruthState.h"
#include "TAxis.h"
#include "TGraph.h"

//=============================================
//======= TFCSEnergyInterpolation =========
//=============================================

TFCSEnergyInterpolationHistogram::TFCSEnergyInterpolationHistogram(
    const char* name, const char* title)
    : TFCSParametrization(name, title)
{
}

auto TFCSEnergyInterpolationHistogram::simulate(
    TFCSSimulationState& simulstate,
    const TFCSTruthState* truth,
    const TFCSExtrapolationState* /*extrapol*/) const -> FCSReturnCode
{
  const double Ekin = truth->Ekin();
  const double Einit = OnlyScaleEnergy() ? simulstate.E() : Ekin;
  // The response is looked up at the true kinetic energy, also when it only
  // scales the energy of the state, as in the other energy interpolations.
  // Outside of the histogram range, the first or last bin is used.
  const int bin = [&] -> int
  {
    if (Ekin < m_hist.GetXaxis()->GetBinLowEdge(1)) {
      return 1;
    }
    if (Ekin > m_hist.GetXaxis()->GetBinUpEdge(m_hist.GetNbinsX())) {
      return m_hist.GetNbinsX();
    }
    return m_hist.GetXaxis()->FindBin(Ekin);
  }();
  const double Emean = m_hist.GetBinContent(bin) * Einit;

  if (OnlyScaleEnergy()) {
    FCS_MSG_DEBUG("set E=" << Emean << " for true Ekin=" << truth->Ekin()
                           << " and E=" << Einit);
  } else {
    FCS_MSG_DEBUG("set E=" << Emean << " for true Ekin=" << truth->Ekin());
  }
  simulstate.set_E(Emean);

  return FCSSuccess;
}

void TFCSEnergyInterpolationHistogram::Print(Option_t* option) const
{
  const auto [shortprint, longprint, optprint] = print_options(option);
  TFCSParametrization::Print(option);

  if (longprint) {
    FCS_MSG_INFO(
        optprint << (OnlyScaleEnergy() ? "  E()*" : "  Ekin()*")
                 << "histNbins=" << m_hist.GetNbinsX() << " "
                 << m_hist.GetXaxis()->GetBinLowEdge(1) << "<=Ekin<="
                 << m_hist.GetXaxis()->GetBinUpEdge(m_hist.GetNbinsX()));
  }
}
