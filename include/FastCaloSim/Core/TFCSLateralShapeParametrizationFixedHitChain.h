// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#ifndef TFCSLateralShapeParametrizationFixedHitChain_h
#define TFCSLateralShapeParametrizationFixedHitChain_h

#include "FastCaloSim/Core/TFCSLateralShapeParametrizationHitChain.h"

class TFCSLateralShapeParametrizationFixedHitChain
    : public TFCSLateralShapeParametrizationHitChain
{
public:
  TFCSLateralShapeParametrizationFixedHitChain(const char* name = nullptr,
                                               const char* title = nullptr);
  TFCSLateralShapeParametrizationFixedHitChain(
      TFCSLateralShapeParametrizationHitBase* hitsim);

protected:
  virtual bool check_all_hits_simulated(
      TFCSLateralShapeParametrizationHitBase::Hit& hit,
      TFCSSimulationState& simulstate,
      const TFCSTruthState* truth,
      const TFCSExtrapolationState* extrapol,
      bool success) const override;

private:
  ClassDefOverride(TFCSLateralShapeParametrizationFixedHitChain,
                   1)  // TFCSLateralShapeParametrizationFixedHitChain
};

#endif
