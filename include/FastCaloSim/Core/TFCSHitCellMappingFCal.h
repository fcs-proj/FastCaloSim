/*
  Copyright (C) 2002-2018 CERN for the benefit of the ATLAS collaboration
*/

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

  virtual FCSReturnCode simulate_hit(
      Hit& hit,
      TFCSSimulationState& simulstate,
      const TFCSTruthState* truth,
      const TFCSExtrapolationState* extrapol) override;

private:
  ClassDefOverride(TFCSHitCellMappingFCal, 1)  // TFCSHitCellMapping
};

#endif
