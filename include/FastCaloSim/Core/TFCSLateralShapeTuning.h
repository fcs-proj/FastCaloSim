// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#ifndef TFCSLateralShapeTuning_h
#define TFCSLateralShapeTuning_h

/*Parent class*/
#include "FastCaloSim/Core/TFCSLateralShapeParametrizationHitBase.h"
/*Truth state include*/
#include "FastCaloSim/Core/TFCSTruthState.h"

/*Pieciewise linear interpolation class include*/
#include "FastCaloSim/Core/TFCSEnergyInterpolationPiecewiseLinear.h"

class TFCSLateralShapeTuning : public TFCSLateralShapeParametrizationHitBase
{
  using interpolationMap =
      std::map<std::string, TFCSEnergyInterpolationPiecewiseLinear*>;

public:
  /// Constructor
  TFCSLateralShapeTuning(const char* name = nullptr,
                         const char* title = nullptr);
  // Destructor
  ~TFCSLateralShapeTuning() override;
  // Simulation call
  auto simulate_hit(Hit& hit,
                    TFCSSimulationState& simulstate,
                    const TFCSTruthState* truth,
                    const TFCSExtrapolationState* extrapol)
      -> FCSReturnCode override;
  // Initialization from model parameter file
  auto initFromModelFile(const std::string& pathToModelParameters,
                         int intMinEta,
                         int intMaxEta) -> FCSReturnCode;
  // Initialization from interpolation map
  auto initFromMap(const interpolationMap&) -> FCSReturnCode;
  // Helper function for strip layer model
  static auto getSeriesScalingFactor(double a0,
                                     double a1,
                                     double a2,
                                     double a3,
                                     double distToShowerCenter) -> double;

private:
  // mapping between the parameter names of the model and its pieciewise linear
  // interpolation object
  interpolationMap m_parameterInterpol;

  ClassDefOverride(TFCSLateralShapeTuning, 1)  // TFCSLateralShapeTuning
};

#endif
