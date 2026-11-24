// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#include <cmath>
#include <cstddef>
#include <iostream>
#include <ranges>
#include <span>
#include <vector>

#include "FastCaloSim/Core/TFCSEnergyInterpolationSpline.h"

#include "FastCaloSim/Core/TFCSExtrapolationState.h"
#include "FastCaloSim/Core/TFCSSimulationState.h"
#include "FastCaloSim/Core/TFCSTruthState.h"

namespace Gaudi::Units
{
constexpr double megaelectronvolt = 1.;
constexpr double kiloelectronvolt = 1.e-3 * megaelectronvolt;
constexpr double keV = kiloelectronvolt;
}  // namespace Gaudi::Units

//=============================================
//======= TFCSEnergyInterpolation =========
//=============================================

TFCSEnergyInterpolationSpline::TFCSEnergyInterpolationSpline(const char* name,
                                                             const char* title)
    : TFCSParametrization(name, title)
{
}

void TFCSEnergyInterpolationSpline::InitFromArrayInLogEkin(Int_t np,
                                                           Double_t* logEkin,
                                                           Double_t* response,
                                                           const char* opt,
                                                           Double_t valbeg,
                                                           Double_t valend)
{
  TSpline3 const initspline(
      GetName(), logEkin, response, np, opt, valbeg, valend);
  m_spline = initspline;
}

void TFCSEnergyInterpolationSpline::InitFromArrayInEkin(Int_t np,
                                                        Double_t* Ekin,
                                                        Double_t* response,
                                                        const char* opt,
                                                        Double_t valbeg,
                                                        Double_t valend)
{
  std::vector<Double_t> logEkin = std::span(Ekin, static_cast<std::size_t>(np))
      | std::views::transform([](Double_t energy) -> Double_t
                              { return std::log(energy); })
      | std::ranges::to<std::vector>();
  InitFromArrayInLogEkin(np, logEkin.data(), response, opt, valbeg, valend);
}

auto TFCSEnergyInterpolationSpline::simulate(
    TFCSSimulationState& simulstate,
    const TFCSTruthState* truth,
    const TFCSExtrapolationState* /*extrapol*/) const -> FCSReturnCode
{
  const double Ekin = truth->Ekin();
  const double Einit = OnlyScaleEnergy() ? simulstate.E() : Ekin;
  // catch very small values of Ekin (use 1 keV here) and fix the spline lookup
  // to the 1keV value
  const double logEkin =
      (Ekin > Gaudi::Units::keV ? std::log(Ekin) : std::log(Gaudi::Units::keV));
  const double Emean = [&] -> double
  {
    if (logEkin < m_spline.GetXmin()) {
      return m_spline.Eval(m_spline.GetXmin()) * Einit;
    }
    if (logEkin > m_spline.GetXmax()) {
      return (m_spline.Eval(m_spline.GetXmax())
              + m_spline.Derivative(m_spline.GetXmax())
                  * (logEkin - m_spline.GetXmax()))
          * Einit;
    }
    return m_spline.Eval(logEkin) * Einit;
  }();

  if (OnlyScaleEnergy()) {
    FCS_MSG_DEBUG("set E=" << Emean << " for true Ekin=" << truth->Ekin()
                           << " and E=" << Einit);
  } else {
    FCS_MSG_DEBUG("set E=" << Emean << " for true Ekin=" << truth->Ekin());
  }
  simulstate.set_E(Emean);

  return FCSSuccess;
}

void TFCSEnergyInterpolationSpline::Print(Option_t* option) const
{
  const auto [shortprint, longprint, optprint] = print_options(option);
  TFCSParametrization::Print(option);

  if (longprint) {
    FCS_MSG_INFO(optprint << (OnlyScaleEnergy() ? "  E()*" : "  Ekin()*")
                          << "Spline N=" << m_spline.GetNp() << " "
                          << m_spline.GetXmin()
                          << "<=log(Ekin)<=" << m_spline.GetXmax() << " "
                          << std::exp(m_spline.GetXmin())
                          << "<=Ekin<=" << std::exp(m_spline.GetXmax()));
  }
}
