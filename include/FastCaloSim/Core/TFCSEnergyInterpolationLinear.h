// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#ifndef ISF_FASTCALOSIMEVENT_TFCSEnergyInterpolationLinear_h
#define ISF_FASTCALOSIMEVENT_TFCSEnergyInterpolationLinear_h

#include "FastCaloSim/Core/TFCSParametrization.h"

class TFCSEnergyInterpolationLinear : public TFCSParametrization
{
public:
  TFCSEnergyInterpolationLinear(const char* name = nullptr,
                                const char* title = nullptr);

  auto is_match_Ekin_bin(int /*Ekin_bin*/) const -> bool override
  {
    return true;
  };
  auto is_match_calosample(int /*calosample*/) const -> bool override
  {
    return true;
  };

  void set_slope(float slope) { m_slope = slope; };
  void set_offset(float offset) { m_offset = offset; };

  // Initialize simulstate with the mean reconstructed energy in the calorimater
  // expected from the true kinetic energy
  auto simulate(TFCSSimulationState& simulstate,
                const TFCSTruthState* truth,
                const TFCSExtrapolationState* extrapol) const
      -> FCSReturnCode override;

  void Print(Option_t* option = "") const override;

private:
  float m_slope;
  float m_offset;

  ClassDefOverride(TFCSEnergyInterpolationLinear,
                   1)  // TFCSEnergyInterpolationLinear
};

#endif
