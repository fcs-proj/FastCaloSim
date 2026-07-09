// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#pragma once

#include "FastCaloSim/Core/TFCS2DFunctionHistogram.h"
#include "FastCaloSim/Core/TFCSLateralShapeParametrizationHitBase.h"
#include "FastCaloSim/Core/TFCSTruthState.h"

class TH2;

class TFCSFlatLateralShapeParametrization
    : public TFCSLateralShapeParametrizationHitBase
{
public:
  TFCSFlatLateralShapeParametrization(const char* name = nullptr,
                                      const char* title = nullptr);
  ~TFCSFlatLateralShapeParametrization() override;

  /// set the integral of the histogram to the desired number of hits
  void set_number_of_hits(float nhits);

  auto get_number_of_expected_hits() const -> float;

  /// default for this class is to simulate poisson(integral histogram) hits
  auto get_number_of_hits(TFCSSimulationState& simulstate,
                          const TFCSTruthState* truth,
                          const TFCSExtrapolationState* extrapol) const
      -> int override;

  /// set the radius in which hits should be generated
  void set_dR(float _dR);

  auto dR() const -> float;

  /// set the radius in which hits should be generated
  void set_scale(float _scale);

  auto scale() const -> float;

  /// simulated one hit position with weight that should be put into simulstate
  /// sometime later all hit weights should be resacled such that their final
  /// sum is simulstate->E(sample) someone also needs to map all hits into cells
  auto simulate_hit(Hit& hit,
                    TFCSSimulationState& simulstate,
                    const TFCSTruthState* truth,
                    const TFCSExtrapolationState* extrapol)
      -> FCSReturnCode override;

  void Print(Option_t* option = "") const override;

protected:
  /// Simulate hits flat in radius dR
  float m_dR {};
  float m_nhits {};
  float m_scale {};

private:
  ClassDefOverride(TFCSFlatLateralShapeParametrization,
                   1)  // TFCSFlatLateralShapeParametrization
};

inline auto TFCSFlatLateralShapeParametrization::get_number_of_expected_hits()
    const -> float
{
  return m_nhits;
}

inline auto TFCSFlatLateralShapeParametrization::dR() const -> float
{
  return m_dR;
}

inline auto TFCSFlatLateralShapeParametrization::scale() const -> float
{
  return m_scale;
}
