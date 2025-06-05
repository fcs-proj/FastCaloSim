// Copyright (c) 2024 CERN for the benefit of the FastCaloSim project

#ifndef TFCSHistoLateralShapeWeightHitAndMiss_h
#define TFCSHistoLateralShapeWeightHitAndMiss_h

#include "FastCaloSim/Core/TFCSHistoLateralShapeWeight.h"

class TH1;
namespace FastCaloSim::Core
{
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

}  // namespace FastCaloSim::Core

#endif
