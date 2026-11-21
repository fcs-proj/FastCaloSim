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
      const TFCSHistoLateralShapeGausLogWeight&) = delete;
  TFCSHistoLateralShapeGausLogWeight(TFCSHistoLateralShapeGausLogWeight&&) =
      delete;
  auto operator=(const TFCSHistoLateralShapeGausLogWeight&)
      -> TFCSHistoLateralShapeGausLogWeight& = delete;
  auto operator=(TFCSHistoLateralShapeGausLogWeight&&)
      -> TFCSHistoLateralShapeGausLogWeight& = delete;

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
