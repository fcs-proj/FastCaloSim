// Copyright (c) 2024 CERN for the benefit of the FastCaloSim project

#ifndef TFCSCenterPositionCalculation_h
#define TFCSCenterPositionCalculation_h

#include "FastCaloSim/Core/TFCSLateralShapeParametrizationHitBase.h"

class TFCSCenterPositionCalculation
    : public TFCSLateralShapeParametrizationHitBase
{
public:
  TFCSCenterPositionCalculation(const char* name = nullptr,
                                const char* title = nullptr);

  /// Used to decorate Hit with extrap center positions
  virtual FCSReturnCode simulate_hit(
      Hit& hit,
      TFCSSimulationState& simulstate,
      const TFCSTruthState* truth,
      const TFCSExtrapolationState* extrapol) override;
  inline void setExtrapWeight(const float weight) { m_extrapWeight = weight; }
  inline float getExtrapWeight() { return m_extrapWeight; }
  void Print(Option_t* option = "") const override;

private:
  float m_extrapWeight;
  ClassDefOverride(TFCSCenterPositionCalculation,
                   1)  // TFCSCenterPositionCalculation
};

#endif
