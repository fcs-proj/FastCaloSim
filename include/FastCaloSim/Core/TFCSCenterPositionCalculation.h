// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#pragma once

#include "FastCaloSim/Core/TFCSLateralShapeParametrizationHitBase.h"

class TFCSCenterPositionCalculation
    : public TFCSLateralShapeParametrizationHitBase
{
public:
  explicit TFCSCenterPositionCalculation(const char* name = nullptr,
                                         const char* title = nullptr);

  /// Used to decorate Hit with extrap center positions
  auto simulate_hit(Hit& hit,
                    TFCSSimulationState& simulstate,
                    const TFCSTruthState* truth,
                    const TFCSExtrapolationState* extrapol)
      -> FCSReturnCode override;
  void setExtrapWeight(const float weight) { m_extrapWeight = weight; }
  auto getExtrapWeight() const -> float { return m_extrapWeight; }
  void Print(Option_t* option = "") const override;

private:
  float m_extrapWeight {0.5};
  ClassDefOverride(TFCSCenterPositionCalculation,
                   1)  // TFCSCenterPositionCalculation
};
