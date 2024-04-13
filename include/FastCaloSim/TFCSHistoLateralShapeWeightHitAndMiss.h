/*
  Copyright (C) 2002-2019 CERN for the benefit of the ATLAS collaboration
*/

#ifndef TFCSHistoLateralShapeWeightHitAndMiss_h
#define TFCSHistoLateralShapeWeightHitAndMiss_h

#include "FastCaloSim/TFCSHistoLateralShapeWeight.h"

class TH1;

class TFCSHistoLateralShapeWeightHitAndMiss : public TFCSHistoLateralShapeWeight
{
public:
  TFCSHistoLateralShapeWeightHitAndMiss(const char* name = nullptr,
                                        const char* title = nullptr);
  virtual ~TFCSHistoLateralShapeWeightHitAndMiss();

  /// weight the energy of one hit in order to generate fluctuations. If the hit
  /// energy is 0, discard the hit
  virtual FCSReturnCode simulate_hit(
      Hit& hit,
      TFCSSimulationState& simulstate,
      const TFCSTruthState* truth,
      const TFCSExtrapolationState* extrapol) override;

private:
  ClassDefOverride(TFCSHistoLateralShapeWeightHitAndMiss,
                   1)  // TFCSHistoLateralShapeWeightHitAndMiss
};

#endif
