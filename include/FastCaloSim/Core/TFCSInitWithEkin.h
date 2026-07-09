// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#pragma once

#include "FastCaloSim/Core/TFCSParametrization.h"

class TFCSInitWithEkin : public TFCSParametrization
{
public:
  TFCSInitWithEkin(const char* name = nullptr, const char* title = nullptr);

  auto is_match_Ekin_bin(int /*Ekin_bin*/) const -> bool override
  {
    return true;
  };
  auto is_match_calosample(int /*calosample*/) const -> bool override
  {
    return true;
  };
  auto is_match_all_Ekin_bin() const -> bool override { return true; };
  auto is_match_all_calosample() const -> bool override { return true; };

  // Initialize simulstate with the kinetic energy Ekin from truth
  auto simulate(TFCSSimulationState& simulstate,
                const TFCSTruthState* truth,
                const TFCSExtrapolationState* extrapol) const
      -> FCSReturnCode override;

private:
  ClassDefOverride(TFCSInitWithEkin, 1)  // TFCSInitWithEkin
};
