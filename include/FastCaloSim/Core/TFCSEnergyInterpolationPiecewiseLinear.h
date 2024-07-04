/*
  Copyright (C) 2002-2020 CERN for the benefit of the ATLAS collaboration
*/

#ifndef ISF_FASTCALOSIMEVENT_TFCSEnergyInterpolationPiecewiseLinear_h
#define ISF_FASTCALOSIMEVENT_TFCSEnergyInterpolationPiecewiseLinear_h

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

  bool OnlyScaleEnergy() const { return TestBit(kOnlyScaleEnergy); };
  void set_OnlyScaleEnergy() { SetBit(kOnlyScaleEnergy); };
  void reset_OnlyScaleEnergy() { ResetBit(kOnlyScaleEnergy); };

  virtual bool is_match_Ekin_bin(int /*Ekin_bin*/) const override
  {
    return true;
  };
  virtual bool is_match_calosample(int /*calosample*/) const override
  {
    return true;
  };

  void InitFromArrayInLogEkin(Int_t np,
                              const Double_t logEkin[],
                              const Double_t response[]);
  void InitFromArrayInEkin(Int_t np,
                           const Double_t Ekin[],
                           const Double_t response[]);

  virtual FCSReturnCode simulate(
      TFCSSimulationState& simulstate,
      const TFCSTruthState* truth,
      const TFCSExtrapolationState* extrapol) const override;

  double evaluate(const double& Ekin) const;

  void Print(Option_t* option = "") const override;

private:
  ROOT::Math::Interpolator m_linInterpol;  //! Do not persistify

  std::vector<double> m_logEkin;
  std::vector<double> m_response;
  std::pair<double, double> m_MinMaxlogEkin;

  ClassDefOverride(TFCSEnergyInterpolationPiecewiseLinear,
                   4)  // TFCSEnergyInterpolationPiecewiseLinear
};

#endif
