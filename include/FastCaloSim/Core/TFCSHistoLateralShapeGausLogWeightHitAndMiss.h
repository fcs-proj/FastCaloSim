// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#pragma once

#include "FastCaloSim/Core/TFCSHistoLateralShapeWeight.h"

class TH1;

class TFCSHistoLateralShapeGausLogWeightHitAndMiss
    : public TFCSHistoLateralShapeWeight
{
public:
  explicit TFCSHistoLateralShapeGausLogWeightHitAndMiss(
      const char* name = nullptr, const char* title = nullptr);
  ~TFCSHistoLateralShapeGausLogWeightHitAndMiss() override;
  TFCSHistoLateralShapeGausLogWeightHitAndMiss(
      const TFCSHistoLateralShapeGausLogWeightHitAndMiss&) = default;
  TFCSHistoLateralShapeGausLogWeightHitAndMiss(
      TFCSHistoLateralShapeGausLogWeightHitAndMiss&&) = default;
  auto operator=(const TFCSHistoLateralShapeGausLogWeightHitAndMiss&)
      -> TFCSHistoLateralShapeGausLogWeightHitAndMiss& = default;
  auto operator=(TFCSHistoLateralShapeGausLogWeightHitAndMiss&&)
      -> TFCSHistoLateralShapeGausLogWeightHitAndMiss& = default;

  /// weight the energy of one hit in order to generate fluctuations. If the hit
  /// energy is 0, discard the hit
  auto simulate_hit(Hit& hit,
                    TFCSSimulationState& simulstate,
                    const TFCSTruthState* truth,
                    const TFCSExtrapolationState* extrapol)
      -> FCSReturnCode override;

private:
  ClassDefOverride(TFCSHistoLateralShapeGausLogWeightHitAndMiss,
                   1)  // TFCSHistoLateralShapeGausLogWeightHitAndMiss
};
