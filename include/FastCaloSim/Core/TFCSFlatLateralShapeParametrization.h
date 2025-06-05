// Copyright (c) 2024 CERN for the benefit of the FastCaloSim project

#ifndef TFCSFlatLateralShapeParametrization_h
#define TFCSFlatLateralShapeParametrization_h

#include "FastCaloSim/Core/TFCS2DFunctionHistogram.h"
#include "FastCaloSim/Core/TFCSLateralShapeParametrizationHitBase.h"
#include "FastCaloSim/Core/TFCSTruthState.h"

class TH2;
namespace FastCaloSim::Core
{
class TFCSFlatLateralShapeParametrization
    : public TFCSLateralShapeParametrizationHitBase
{
public:
  TFCSFlatLateralShapeParametrization(const char* name = nullptr,
                                      const char* title = nullptr);
  virtual ~TFCSFlatLateralShapeParametrization();

  /// set the integral of the histogram to the desired number of hits
  void set_number_of_hits(float nhits);

  float get_number_of_expected_hits() const;

  /// default for this class is to simulate poisson(integral histogram) hits
  virtual int get_number_of_hits(
      TFCSSimulationState& simulstate,
      const TFCSTruthState* truth,
      const TFCSExtrapolationState* extrapol) const override;

  /// set the radius in which hits should be generated
  void set_dR(float _dR);

  float dR() const;

  /// set the radius in which hits should be generated
  void set_scale(float _scale);

  float scale() const;

  /// simulated one hit position with weight that should be put into simulstate
  /// sometime later all hit weights should be resacled such that their final
  /// sum is simulstate->E(sample) someone also needs to map all hits into cells
  virtual FCSReturnCode simulate_hit(
      Hit& hit,
      TFCSSimulationState& simulstate,
      const TFCSTruthState* truth,
      const TFCSExtrapolationState* extrapol) override;

  virtual void Print(Option_t* option = "") const override;

protected:
  /// Simulate hits flat in radius dR
  float m_dR;
  float m_nhits;
  float m_scale;

private:
  ClassDefOverride(TFCSFlatLateralShapeParametrization,
                   1)  // TFCSFlatLateralShapeParametrization
};

inline float TFCSFlatLateralShapeParametrization::get_number_of_expected_hits()
    const
{
  return m_nhits;
}

inline float TFCSFlatLateralShapeParametrization::dR() const
{
  return m_dR;
}

inline float TFCSFlatLateralShapeParametrization::scale() const
{
  return m_scale;
}
}  // namespace FastCaloSim::Core
#endif
