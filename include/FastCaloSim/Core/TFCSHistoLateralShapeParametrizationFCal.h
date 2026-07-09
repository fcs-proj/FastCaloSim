// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#pragma once

#include "FastCaloSim/Core/TFCS2DFunctionHistogram.h"
#include "FastCaloSim/Core/TFCSHistoLateralShapeParametrization.h"
#include "FastCaloSim/Core/TFCSTruthState.h"

class TH2;

class TFCSHistoLateralShapeParametrizationFCal
    : public TFCSHistoLateralShapeParametrization
{
public:
  TFCSHistoLateralShapeParametrizationFCal(const char* name = nullptr,
                                           const char* title = nullptr);
  ~TFCSHistoLateralShapeParametrizationFCal() override;

  auto simulate_hit(Hit& hit,
                    TFCSSimulationState& simulstate,
                    const TFCSTruthState* truth,
                    const TFCSExtrapolationState* extrapol)
      -> FCSReturnCode override;

private:
  ClassDefOverride(TFCSHistoLateralShapeParametrizationFCal,
                   1)  // TFCSHistoLateralShapeParametrizationFCal
};
