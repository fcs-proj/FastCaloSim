// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#pragma once

#include "FastCaloSim/Core/TFCSHistoLateralShapeWeight.h"

class TH1;

class TFCSHistoLateralShapeGausLogWeight : public TFCSHistoLateralShapeWeight
{
public:
  explicit TFCSHistoLateralShapeGausLogWeight(const char* name = nullptr,
                                              const char* title = nullptr);
  ~TFCSHistoLateralShapeGausLogWeight() override;
  TFCSHistoLateralShapeGausLogWeight(
      const TFCSHistoLateralShapeGausLogWeight&) = default;
  TFCSHistoLateralShapeGausLogWeight(TFCSHistoLateralShapeGausLogWeight&&) =
      default;
  auto operator=(const TFCSHistoLateralShapeGausLogWeight&)
      -> TFCSHistoLateralShapeGausLogWeight& = default;
  auto operator=(TFCSHistoLateralShapeGausLogWeight&&)
      -> TFCSHistoLateralShapeGausLogWeight& = default;

  /// weight the energy of one hit in order to generate fluctuations
  auto simulate_hit(Hit& hit,
                    TFCSSimulationState& simulstate,
                    const TFCSTruthState* truth,
                    const TFCSExtrapolationState* extrapol)
      -> FCSReturnCode override;

protected:
ClassDefOverride(TFCSHistoLateralShapeGausLogWeight,
                 1)  // TFCSHistoLateralShapeGausLogWeight

    private:
};
