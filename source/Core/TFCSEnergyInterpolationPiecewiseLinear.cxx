// Copyright (c) 2024 CERN for the benefit of the FastCaloSim project

#include "FastCaloSim/Core/TFCSEnergyInterpolationPiecewiseLinear.h"

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

//========================================================
//======= TFCSEnergyInterpolationPiecewiseLinear =========
//========================================================

TFCSEnergyInterpolationPiecewiseLinear::TFCSEnergyInterpolationPiecewiseLinear(
    const char* name, const char* title)
    : TFCSParametrization(name, title)
    , m_linInterpol(0, ROOT::Math::Interpolation::kLINEAR)
{
}

void TFCSEnergyInterpolationPiecewiseLinear::InitFromArrayInLogEkin(
    Int_t np, const Double_t logEkin[], const Double_t response[])
{
  // save logEkin and response as std::vector class members
  // this is required for using the custom streamer
  m_logEkin.assign(logEkin, logEkin + np);
  m_response.assign(response, response + np);
  m_linInterpol.SetData(m_logEkin, m_response);

  auto min_max = std::minmax_element(m_logEkin.begin(), m_logEkin.end());
  m_MinMaxlogEkin = std::make_pair(*min_max.first, *min_max.second);
}

void TFCSEnergyInterpolationPiecewiseLinear::InitFromArrayInEkin(
    Int_t np, const Double_t Ekin[], const Double_t response[])
{
  std::vector<Double_t> logEkin(np);
  for (int i = 0; i < np; i++)
    logEkin[i] = TMath::Log(Ekin[i]);
  InitFromArrayInLogEkin(np, logEkin.data(), response);
}

FCSReturnCode TFCSEnergyInterpolationPiecewiseLinear::simulate(
    TFCSSimulationState& simulstate,
    const TFCSTruthState* truth,
    const TFCSExtrapolationState*) const
{
  const float Ekin = truth->Ekin();
  const float Einit = OnlyScaleEnergy() ? simulstate.E() : Ekin;

  // catch very small values of Ekin (use 1 keV here) and fix the interpolation
  // lookup to the 1keV value
  const float logEkin = Ekin > Gaudi::Units::keV
      ? TMath::Log(Ekin)
      : TMath::Log(Gaudi::Units::keV);

  float Emean;
  if (logEkin < m_MinMaxlogEkin.first) {
    Emean = m_linInterpol.Eval(m_MinMaxlogEkin.first) * Einit;
  } else if (logEkin > m_MinMaxlogEkin.second) {
    Emean = (m_linInterpol.Eval(m_MinMaxlogEkin.second)
             + m_linInterpol.Deriv(m_MinMaxlogEkin.second)
                 * (logEkin - m_MinMaxlogEkin.second))
        * Einit;
  } else {
    Emean = m_linInterpol.Eval(logEkin) * Einit;
  }

  if (OnlyScaleEnergy())
    FCS_MSG_DEBUG("set E=" << Emean << " for true Ekin=" << Ekin
                           << " and E=" << Einit);
  else
    FCS_MSG_DEBUG("set E=" << Emean << " for true Ekin=" << Ekin);

  // set mean energy of simulstate
  simulstate.set_E(Emean);

  return FCSSuccess;
}

double TFCSEnergyInterpolationPiecewiseLinear::evaluate(
    const double& Ekin) const
{
  // returns simple evaluation of the interpolation
  // if the lookup is below the minimum interpolation value, will return minimum
  // evaluation if the lookup is above the maximum interpolation value, will
  // return maximum evaluation this means no extrapolation beyond the
  // interpolation will be performed

  // catch very small values of Ekin (use 1 keV here) and fix the interpolation
  // lookup to the 1keV value
  const float logEkin = Ekin > Gaudi::Units::keV
      ? TMath::Log(Ekin)
      : TMath::Log(Gaudi::Units::keV);

  if (logEkin < m_MinMaxlogEkin.first)
    return m_linInterpol.Eval(m_MinMaxlogEkin.first);
  else if (logEkin > m_MinMaxlogEkin.second)
    return m_linInterpol.Eval(m_MinMaxlogEkin.second);
  else
    return m_linInterpol.Eval(logEkin);
}

void TFCSEnergyInterpolationPiecewiseLinear::Print(Option_t* option) const
{
  TString opt(option);
  bool shortprint = opt.Index("short") >= 0;
  bool longprint =
      msgLvl(FCS_MSG::DEBUG) || (msgLvl(FCS_MSG::INFO) && !shortprint);
  TString optprint = opt;
  optprint.ReplaceAll("short", "");
  TFCSParametrization::Print(option);

  if (longprint)
    FCS_MSG_INFO(optprint << (OnlyScaleEnergy() ? "  E()*" : "  Ekin()*")
                          << "linInterpol N=" << m_logEkin.size() << " "
                          << m_MinMaxlogEkin.first
                          << "<=log(Ekin)<=" << m_MinMaxlogEkin.second << " "
                          << TMath::Exp(m_MinMaxlogEkin.first)
                          << "<=Ekin<=" << TMath::Exp(m_MinMaxlogEkin.second));
}

void TFCSEnergyInterpolationPiecewiseLinear::Streamer(TBuffer& R__b)
{
  // Stream an object of class TFCSEnergyInterpolationPiecewiseLinear
  if (R__b.IsReading()) {
    // read the class buffer
    R__b.ReadClassBuffer(TFCSEnergyInterpolationPiecewiseLinear::Class(), this);
    // initialize interpolation from saved class members
    InitFromArrayInLogEkin(
        m_logEkin.size(), m_logEkin.data(), m_response.data());
  } else {
    R__b.WriteClassBuffer(TFCSEnergyInterpolationPiecewiseLinear::Class(),
                          this);
  }
}
