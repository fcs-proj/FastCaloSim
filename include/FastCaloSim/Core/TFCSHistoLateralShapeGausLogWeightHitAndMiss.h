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
      const TFCSHistoLateralShapeGausLogWeightHitAndMiss&) = delete;
  TFCSHistoLateralShapeGausLogWeightHitAndMiss(
      TFCSHistoLateralShapeGausLogWeightHitAndMiss&&) = delete;
  auto operator=(const TFCSHistoLateralShapeGausLogWeightHitAndMiss&)
      -> TFCSHistoLateralShapeGausLogWeightHitAndMiss& = delete;
  auto operator=(TFCSHistoLateralShapeGausLogWeightHitAndMiss&&)
      -> TFCSHistoLateralShapeGausLogWeightHitAndMiss& = delete;

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
