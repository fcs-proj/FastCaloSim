/*
  Copyright (C) 2002-2018 CERN for the benefit of the ATLAS collaboration
*/

#ifndef TFCSHistoLateralShapeParametrizationFCal_h
#define TFCSHistoLateralShapeParametrizationFCal_h

#include "FastCaloSim/TFCSHistoLateralShapeParametrization.h"
#include "FastCaloSim/TFCS2DFunctionHistogram.h"
#include "FastCaloSim/TFCSTruthState.h"

class TH2;

class TFCSHistoLateralShapeParametrizationFCal
    : public TFCSHistoLateralShapeParametrization {
public:
  TFCSHistoLateralShapeParametrizationFCal(const char *name = nullptr,
                                           const char *title = nullptr);
  ~TFCSHistoLateralShapeParametrizationFCal();

  virtual FCSReturnCode
  simulate_hit(Hit &hit, TFCSSimulationState &simulstate,
               const TFCSTruthState *truth,
               const TFCSExtrapolationState *extrapol) override;

private:
  ClassDefOverride(TFCSHistoLateralShapeParametrizationFCal,
                   1) // TFCSHistoLateralShapeParametrizationFCal
};

#endif
