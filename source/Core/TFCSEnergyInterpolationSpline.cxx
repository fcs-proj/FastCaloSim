/*
  Copyright (C) 2002-2020 CERN for the benefit of the ATLAS collaboration
*/

#include <iostream>
#include <vector>

#include "FastCaloSim/Core/TFCSEnergyInterpolationSpline.h"

#include "FastCaloSim/Core/TFCSExtrapolationState.h"
#include "FastCaloSim/Core/TFCSSimulationState.h"
#include "FastCaloSim/Core/TFCSTruthState.h"

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

TFCSEnergyInterpolationSpline::TFCSEnergyInterpolationSpline(const char* name,
                                                             const char* title)
    : TFCSParametrization(name, title)
{
}

void TFCSEnergyInterpolationSpline::InitFromArrayInLogEkin(Int_t np,
                                                           Double_t logEkin[],
                                                           Double_t response[],
                                                           const char* opt,
                                                           Double_t valbeg,
                                                           Double_t valend)
{
  TSpline3 initspline(GetName(), logEkin, response, np, opt, valbeg, valend);
  m_spline = initspline;
}

void TFCSEnergyInterpolationSpline::InitFromArrayInEkin(Int_t np,
                                                        Double_t Ekin[],
                                                        Double_t response[],
                                                        const char* opt,
                                                        Double_t valbeg,
                                                        Double_t valend)
{
  std::vector<Double_t> logEkin(np);
  for (int i = 0; i < np; ++i)
    logEkin[i] = TMath::Log(Ekin[i]);
  InitFromArrayInLogEkin(np, logEkin.data(), response, opt, valbeg, valend);
}

FCSReturnCode TFCSEnergyInterpolationSpline::simulate(
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
  // catch very small values of Ekin (use 1 keV here) and fix the spline lookup
  // to the 1keV value
  const float logEkin =
      (Ekin > Gaudi::Units::keV ? TMath::Log(Ekin)
                                : TMath::Log(Gaudi::Units::keV));
  if (logEkin < m_spline.GetXmin()) {
    Emean = m_spline.Eval(m_spline.GetXmin()) * Einit;
  } else {
    if (logEkin > m_spline.GetXmax()) {
      Emean = (m_spline.Eval(m_spline.GetXmax())
               + m_spline.Derivative(m_spline.GetXmax())
                   * (logEkin - m_spline.GetXmax()))
          * Einit;
    } else {
      Emean = m_spline.Eval(logEkin) * Einit;
    }
  }

  if (OnlyScaleEnergy()) {
    ATH_MSG_DEBUG("set E=" << Emean << " for true Ekin=" << truth->Ekin()
                           << " and E=" << Einit);
  } else {
    ATH_MSG_DEBUG("set E=" << Emean << " for true Ekin=" << truth->Ekin());
  }
  simulstate.set_E(Emean);

  return FCSSuccess;
}

void TFCSEnergyInterpolationSpline::Print(Option_t* option) const
{
  TString opt(option);
  bool shortprint = opt.Index("short") >= 0;
  bool longprint = msgLvl(MSG::DEBUG) || (msgLvl(MSG::INFO) && !shortprint);
  TString optprint = opt;
  optprint.ReplaceAll("short", "");
  TFCSParametrization::Print(option);

  if (longprint)
    ATH_MSG_INFO(optprint << (OnlyScaleEnergy() ? "  E()*" : "  Ekin()*")
                          << "Spline N=" << m_spline.GetNp() << " "
                          << m_spline.GetXmin()
                          << "<=log(Ekin)<=" << m_spline.GetXmax() << " "
                          << TMath::Exp(m_spline.GetXmin())
                          << "<=Ekin<=" << TMath::Exp(m_spline.GetXmax()));
}
