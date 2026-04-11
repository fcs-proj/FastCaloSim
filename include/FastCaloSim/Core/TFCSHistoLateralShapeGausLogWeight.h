// Copyright (c) 2024 CERN for the benefit of the FastCaloSim project

#ifndef TFCSHistoLateralShapeGausLogWeight_h
#define TFCSHistoLateralShapeGausLogWeight_h

#include "FastCaloSim/Core/TFCSHistoLateralShapeWeight.h"

class TH1;
namespace FastCaloSim::Core
{
class TFCSHistoLateralShapeGausLogWeight : public TFCSHistoLateralShapeWeight
{
public:
  TFCSHistoLateralShapeGausLogWeight(const char* name = nullptr,
                                     const char* title = nullptr);
  virtual ~TFCSHistoLateralShapeGausLogWeight();

  /// weight the energy of one hit in order to generate fluctuations
  virtual FCSReturnCode simulate_hit(
      Hit& hit,
      TFCSSimulationState& simulstate,
      const TFCSTruthState* truth,
      const TFCSExtrapolationState* extrapol) override;

protected:
  ClassDefOverride(TFCSHistoLateralShapeGausLogWeight,
                   1)  // TFCSHistoLateralShapeGausLogWeight

      private:
};
}  // namespace FastCaloSim::Core
#endif
