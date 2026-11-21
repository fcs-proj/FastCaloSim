// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#pragma once

#include "FastCaloSim/Core/TFCSHistoLateralShapeWeight.h"

class TH1;

class TFCSHistoLateralShapeWeightHitAndMiss : public TFCSHistoLateralShapeWeight
{
public:
  explicit TFCSHistoLateralShapeWeightHitAndMiss(const char* name = nullptr,
                                                 const char* title = nullptr);
  ~TFCSHistoLateralShapeWeightHitAndMiss() override;
  TFCSHistoLateralShapeWeightHitAndMiss(
      const TFCSHistoLateralShapeWeightHitAndMiss&) = default;
  TFCSHistoLateralShapeWeightHitAndMiss(
      TFCSHistoLateralShapeWeightHitAndMiss&&) = default;
  auto operator=(const TFCSHistoLateralShapeWeightHitAndMiss&)
      -> TFCSHistoLateralShapeWeightHitAndMiss& = default;
  auto operator=(TFCSHistoLateralShapeWeightHitAndMiss&&)
      -> TFCSHistoLateralShapeWeightHitAndMiss& = default;

  /// weight the energy of one hit in order to generate fluctuations. If the hit
  /// energy is 0, discard the hit
  auto simulate_hit(Hit& hit,
                    TFCSSimulationState& simulstate,
                    const TFCSTruthState* truth,
                    const TFCSExtrapolationState* extrapol)
      -> FCSReturnCode override;

private:
  ClassDefOverride(TFCSHistoLateralShapeWeightHitAndMiss,
                   1)  // TFCSHistoLateralShapeWeightHitAndMiss
};
