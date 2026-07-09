// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#pragma once

#include "FastCaloSim/Core/TFCSParametrization.h"

class TFCSInvisibleParametrization : public TFCSParametrization
{
public:
  TFCSInvisibleParametrization(const char* name = nullptr,
                               const char* title = nullptr)
      : TFCSParametrization(name, title) {};

  auto is_match_Ekin_bin(int /*Ekin_bin*/) const -> bool override
  {
    return true;
  };
  auto is_match_calosample(int /*calosample*/) const -> bool override
  {
    return true;
  };

  auto simulate(TFCSSimulationState& simulstate,
                const TFCSTruthState* truth,
                const TFCSExtrapolationState* extrapol) const
      -> FCSReturnCode override;

private:
  ClassDefOverride(TFCSInvisibleParametrization,
                   1)  // TFCSInvisibleParametrization
};
