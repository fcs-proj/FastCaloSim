// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#pragma once

// base class include
#include "FastCaloSim/Core/TFCSParametrization.h"
// interpolator include
#include <Math/Interpolator.h>
// TBuffer include required for custom class streamer
#include "TBuffer.h"

class TFCSEnergyInterpolationPiecewiseLinear : public TFCSParametrization
{
public:
  TFCSEnergyInterpolationPiecewiseLinear(const char* name = nullptr,
                                         const char* title = nullptr);
  /// Status bit for energy initialization
  enum FCSEnergyInitializationStatusBits
  {
    kOnlyScaleEnergy =
        BIT(15)  ///< Set this bit in the TObject bit field the simulated energy
                 ///< should only be scaled by the spline
  };

  auto OnlyScaleEnergy() const -> bool { return TestBit(kOnlyScaleEnergy); };
  void set_OnlyScaleEnergy() { SetBit(kOnlyScaleEnergy); };
  void reset_OnlyScaleEnergy() { ResetBit(kOnlyScaleEnergy); };

  auto is_match_Ekin_bin(int /*Ekin_bin*/) const -> bool override
  {
    return true;
  };
  auto is_match_calosample(int /*calosample*/) const -> bool override
  {
    return true;
  };

  void InitFromArrayInLogEkin(Int_t np,
                              const Double_t logEkin[],
                              const Double_t response[]);
  void InitFromArrayInEkin(Int_t np,
                           const Double_t Ekin[],
                           const Double_t response[]);

  auto simulate(TFCSSimulationState& simulstate,
                const TFCSTruthState* truth,
                const TFCSExtrapolationState* extrapol) const
      -> FCSReturnCode override;

  auto evaluate(const double& Ekin) const -> double;

  void Print(Option_t* option = "") const override;

private:
  ROOT::Math::Interpolator m_linInterpol;  //! Do not persistify

  std::vector<double> m_logEkin;
  std::vector<double> m_response;
  std::pair<double, double> m_MinMaxlogEkin;

  ClassDefOverride(TFCSEnergyInterpolationPiecewiseLinear,
                   4)  // TFCSEnergyInterpolationPiecewiseLinear
};
