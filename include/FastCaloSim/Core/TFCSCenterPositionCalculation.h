// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#pragma once

#include "FastCaloSim/Core/TFCSLateralShapeParametrizationHitBase.h"

class TFCSCenterPositionCalculation
    : public TFCSLateralShapeParametrizationHitBase
{
public:
  TFCSCenterPositionCalculation(const char* name = nullptr,
                                const char* title = nullptr);

  /// Used to decorate Hit with extrap center positions
  auto simulate_hit(Hit& hit,
                    TFCSSimulationState& simulstate,
                    const TFCSTruthState* truth,
                    const TFCSExtrapolationState* extrapol)
      -> FCSReturnCode override;
  inline void setExtrapWeight(const float weight) { m_extrapWeight = weight; }
  inline auto getExtrapWeight() -> float { return m_extrapWeight; }
  void Print(Option_t* option = "") const override;

private:
  float m_extrapWeight;
  ClassDefOverride(TFCSCenterPositionCalculation,
                   1)  // TFCSCenterPositionCalculation
};
