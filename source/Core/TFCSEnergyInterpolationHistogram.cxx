// Copyright (c) 2024 CERN for the benefit of the FastCaloSim project

#include <iostream>

#include "FastCaloSim/Core/TFCSEnergyInterpolationHistogram.h"

#include "FastCaloSim/Core/TFCSExtrapolationState.h"
#include "FastCaloSim/Core/TFCSSimulationState.h"
#include "FastCaloSim/Core/TFCSTruthState.h"
#include "TAxis.h"
#include "TGraph.h"

namespace Gaudi
{
namespace Units
{
constexpr double megaelectronvolt = 1.;
constexpr double kiloelectronvolt = 1.e-3 * megaelectronvolt;
constexpr double keV = kiloelectronvolt;
}  // namespace Units
}  // namespace Gaudi

//=============================================
//======= TFCSEnergyInterpolation =========
//=============================================

TFCSEnergyInterpolationHistogram::TFCSEnergyInterpolationHistogram(
    const char* name, const char* title)
    : TFCSParametrization(name, title)
{
}

FCSReturnCode TFCSEnergyInterpolationHistogram::simulate(
    TFCSSimulationState& simulstate,
    const TFCSTruthState* truth,
    const TFCSExtrapolationState*) const
{
  float Emean;
  float Einit;
  const float Ekin = truth->Ekin();
  if (OnlyScaleEnergy())
    Einit = simulstate.E();
  else
    Einit = Ekin;
  if (Einit < m_hist.GetXaxis()->GetBinLowEdge(1)) {
    Emean = m_hist.GetBinContent(1) * Einit;
  } else {
    if (Einit > m_hist.GetXaxis()->GetBinUpEdge(m_hist.GetNbinsX())) {
      Emean = m_hist.GetBinContent(m_hist.GetNbinsX()) * Einit;
    } else {
      Emean = m_hist.GetBinContent(m_hist.GetXaxis()->FindBin(Einit)) * Einit;
    }
  }

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
  TString opt(option);
  bool shortprint = opt.Index("short") >= 0;
  bool longprint =
      msgLvl(FCS_MSG::DEBUG) || (msgLvl(FCS_MSG::INFO) && !shortprint);
  TString optprint = opt;
  optprint.ReplaceAll("short", "");
  TFCSParametrization::Print(option);

  if (longprint)
    FCS_MSG_INFO(
        optprint << (OnlyScaleEnergy() ? "  E()*" : "  Ekin()*")
                 << "histNbins=" << m_hist.GetNbinsX() << " "
                 << m_hist.GetXaxis()->GetBinLowEdge(1) << "<=Ekin<="
                 << m_hist.GetXaxis()->GetBinUpEdge(m_hist.GetNbinsX()));
}
