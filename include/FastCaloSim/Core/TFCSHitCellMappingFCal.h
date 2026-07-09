// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#ifndef TFCSHitCellMappingFCal_h
#define TFCSHitCellMappingFCal_h

#include "FastCaloSim/Core/TFCSHitCellMapping.h"

class CaloGeo;

class TFCSHitCellMappingFCal : public TFCSHitCellMapping
{
public:
  TFCSHitCellMappingFCal(const char* name = nullptr,
                         const char* title = nullptr,
                         CaloGeo* geo = nullptr)
      : TFCSHitCellMapping(name, title, geo)
  {
  }

  auto simulate_hit(Hit& hit,
                    TFCSSimulationState& simulstate,
                    const TFCSTruthState* truth,
                    const TFCSExtrapolationState* extrapol)
      -> FCSReturnCode override;

private:
  ClassDefOverride(TFCSHitCellMappingFCal, 1)  // TFCSHitCellMapping
};

#endif
