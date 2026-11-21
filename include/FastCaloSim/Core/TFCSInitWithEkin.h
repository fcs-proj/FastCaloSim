// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#pragma once

#include "FastCaloSim/Core/TFCSParametrization.h"

class TFCSInitWithEkin : public TFCSParametrization
{
public:
  explicit TFCSInitWithEkin(const char* name = nullptr,
                            const char* title = nullptr);

  auto is_match_Ekin_bin(std::size_t /*Ekin_bin*/) const -> bool override
  {
    return true;
  }
  auto is_match_calosample(std::size_t /*calosample*/) const -> bool override
  {
    return true;
  }
  auto is_match_all_Ekin_bin() const -> bool override { return true; }
  auto is_match_all_calosample() const -> bool override { return true; }

  // Initialize simulstate with the kinetic energy Ekin from truth
  auto simulate(TFCSSimulationState& simulstate,
                const TFCSTruthState* truth,
                const TFCSExtrapolationState* extrapol) const
      -> FCSReturnCode override;

private:
  ClassDefOverride(TFCSInitWithEkin, 1)  // TFCSInitWithEkin
};
