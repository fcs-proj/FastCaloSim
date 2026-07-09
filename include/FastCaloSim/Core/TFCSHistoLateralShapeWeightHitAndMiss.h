// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#ifndef TFCSHistoLateralShapeWeightHitAndMiss_h
#define TFCSHistoLateralShapeWeightHitAndMiss_h

#include "FastCaloSim/Core/TFCSHistoLateralShapeWeight.h"

class TH1;

class TFCSHistoLateralShapeWeightHitAndMiss : public TFCSHistoLateralShapeWeight
{
public:
  TFCSHistoLateralShapeWeightHitAndMiss(const char* name = nullptr,
                                        const char* title = nullptr);
  ~TFCSHistoLateralShapeWeightHitAndMiss() override;

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

#endif
