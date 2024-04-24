/*
  Copyright (C) 2002-2019 CERN for the benefit of the ATLAS collaboration
*/

#ifndef TFCSHistoLateralShapeGausLogWeightHitAndMiss_h
#define TFCSHistoLateralShapeGausLogWeightHitAndMiss_h

#include "FastCaloSim/Core/TFCSHistoLateralShapeWeight.h"

class TH1;

class TFCSHistoLateralShapeGausLogWeightHitAndMiss
    : public TFCSHistoLateralShapeWeight
{
public:
  TFCSHistoLateralShapeGausLogWeightHitAndMiss(const char* name = nullptr,
                                               const char* title = nullptr);
  virtual ~TFCSHistoLateralShapeGausLogWeightHitAndMiss();

  /// weight the energy of one hit in order to generate fluctuations. If the hit
  /// energy is 0, discard the hit
  virtual FCSReturnCode simulate_hit(
      Hit& hit,
      TFCSSimulationState& simulstate,
      const TFCSTruthState* truth,
      const TFCSExtrapolationState* extrapol) override;

private:
  ClassDefOverride(TFCSHistoLateralShapeGausLogWeightHitAndMiss,
                   1)  // TFCSHistoLateralShapeGausLogWeightHitAndMiss
};

#endif
