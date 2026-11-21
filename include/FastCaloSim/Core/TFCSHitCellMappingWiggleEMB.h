// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#pragma once

#include "FastCaloSim/Core/TFCSHitCellMapping.h"

class TFCSHitCellMappingWiggleEMB : public TFCSHitCellMapping
{
public:
  explicit TFCSHitCellMappingWiggleEMB(const char* name = nullptr,
                                       const char* title = nullptr,
                                       CaloGeo* geo = nullptr);

  /// modify one hit position to emulate the LAr accordion shape
  /// and then fills all hits into calorimeter cells
  auto simulate_hit(Hit& hit,
                    TFCSSimulationState& simulstate,
                    const TFCSTruthState* truth,
                    const TFCSExtrapolationState* extrapol)
      -> FCSReturnCode override;

private:
  //** Array for the hit-to-cell assignment accordion structure fix (wiggle)
  //**//
  //** To be moved to the conditions database at some point **//
  // Persistent members: the layout is part of the file format
  // NOLINTBEGIN(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays)
  double m_wiggleLayer1[50];
  double m_wiggleLayer2[50];
  double m_wiggleLayer3[50];
  // NOLINTEND(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays)

  auto doWiggle(double searchRand) -> double;

  ClassDefOverride(TFCSHitCellMappingWiggleEMB,
                   1)  // TFCSHitCellMappingWiggleEMB
};
