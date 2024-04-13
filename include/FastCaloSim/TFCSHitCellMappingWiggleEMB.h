/*
  Copyright (C) 2002-2018 CERN for the benefit of the ATLAS collaboration
*/

#ifndef TFCSHitCellMappingWiggleEMB_h
#define TFCSHitCellMappingWiggleEMB_h

#include "FastCaloSim/TFCSHitCellMapping.h"

class TFCSHitCellMappingWiggleEMB : public TFCSHitCellMapping {
public:
  TFCSHitCellMappingWiggleEMB(const char *name = nullptr,
                              const char *title = nullptr,
                              ICaloGeometry *geo = nullptr);

  /// modify one hit position to emulate the LAr accordeon shape
  /// and then fills all hits into calorimeter cells
  virtual FCSReturnCode
  simulate_hit(Hit &hit, TFCSSimulationState &simulstate,
               const TFCSTruthState *truth,
               const TFCSExtrapolationState *extrapol) override;

private:
  //** Array for the hit-to-cell assignment accordion structure fix (wiggle)
  //**//
  //** To be moved to the conditions database at some point **//
  double m_wiggleLayer1[50];
  double m_wiggleLayer2[50];
  double m_wiggleLayer3[50];

  double doWiggle(double searchRand);

  ClassDefOverride(TFCSHitCellMappingWiggleEMB,
                   1) // TFCSHitCellMappingWiggleEMB
};

#endif
