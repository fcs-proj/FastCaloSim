// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#pragma once

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
  explicit TFCSLateralShapeTuning(const char* name = nullptr,
                                  const char* title = nullptr);
  // Destructor
  ~TFCSLateralShapeTuning() override;
  // Owns its pointers: a copy would delete them twice
  TFCSLateralShapeTuning(const TFCSLateralShapeTuning&) = delete;
  TFCSLateralShapeTuning(TFCSLateralShapeTuning&&) = delete;
  auto operator=(const TFCSLateralShapeTuning&)
      -> TFCSLateralShapeTuning& = delete;
  auto operator=(TFCSLateralShapeTuning&&) -> TFCSLateralShapeTuning& = delete;
  // Simulation call: applies the model whose parameters were loaded
  auto simulate_hit(Hit& hit,
                    TFCSSimulationState& simulstate,
                    const TFCSTruthState* truth,
                    const TFCSExtrapolationState* extrapol)
      -> FCSReturnCode override;
  /// The tuning models
  enum class Model
  {
    /// Scaling of the distance to the shower center in eta by a series in
    /// that distance (parameters a0 to a3), for strip layers
    kSeries,
    /// Scaling of the distances to the shower center in eta and phi (eta_s
    /// and phi_s), for middle layers
    kEtaPhiScaling
  };
  // Initialization from model parameter file: the parameters of the given
  // model for the layer of this object
  auto initFromModelFile(const std::string& pathToModelParameters,
                         int intMinEta,
                         int intMaxEta,
                         Model model) -> FCSReturnCode;
  // Initialization from interpolation map
  auto initFromMap(const interpolationMap& /*parameters*/) -> FCSReturnCode;
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
