// Copyright (c) 2024 CERN for the benefit of the FastCaloSim project

#ifndef TFCSHitCellMappingWiggleEMB_h
#define TFCSHitCellMappingWiggleEMB_h

#include "FastCaloSim/Core/TFCSHitCellMapping.h"

namespace FastCaloSim::Core
{

class TFCSHitCellMappingWiggleEMB : public TFCSHitCellMapping
{
public:
  TFCSHitCellMappingWiggleEMB(const char* name = nullptr,
                              const char* title = nullptr,
                              CaloGeo* geo = nullptr);

  /// modify one hit position to emulate the LAr accordion shape
  /// and then fills all hits into calorimeter cells
  virtual FCSReturnCode simulate_hit(
      Hit& hit,
      TFCSSimulationState& simulstate,
      const TFCSTruthState* truth,
      const TFCSExtrapolationState* extrapol) override;

private:
  //** Array for the hit-to-cell assignment accordion structure fix (wiggle)
  //**//
  //** To be moved to the conditions database at some point **//
  double m_wiggleLayer1[50];
  double m_wiggleLayer2[50];
  double m_wiggleLayer3[50];

  double doWiggle(double searchRand);

  ClassDefOverride(TFCSHitCellMappingWiggleEMB,
                   1)  // TFCSHitCellMappingWiggleEMB
};
}  // namespace FastCaloSim::Core
#endif
