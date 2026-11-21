// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#pragma once

#include "FastCaloSim/Core/TFCSLateralShapeParametrizationHitChain.h"

class TFCSLateralShapeParametrizationFluctChain
    : public TFCSLateralShapeParametrizationHitChain
{
public:
  explicit TFCSLateralShapeParametrizationFluctChain(
      const char* name = nullptr, const char* title = nullptr, float RMS = 1.0);
  explicit TFCSLateralShapeParametrizationFluctChain(
      TFCSLateralShapeParametrizationHitBase* hitsim);

  // set and get the amount of hit energy fluctation around E/n for n hits
  void set_hit_RMS(float RMS);
  auto get_hit_RMS() const -> float { return m_RMS; }

  /// Get hit energy from layer energy and number of hits
  auto get_E_hit(TFCSSimulationState& simulstate,
                 const TFCSTruthState* truth,
                 const TFCSExtrapolationState* extrapol) const
      -> double override;

  auto simulate(TFCSSimulationState& simulstate,
                const TFCSTruthState* truth,
                const TFCSExtrapolationState* extrapol) const
      -> FCSReturnCode override;

  void Print(Option_t* option) const override;

private:
  float m_RMS {};

  ClassDefOverride(TFCSLateralShapeParametrizationFluctChain,
                   1)  // TFCSLateralShapeParametrizationFluctChain
};
