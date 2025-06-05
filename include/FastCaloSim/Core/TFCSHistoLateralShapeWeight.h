// Copyright (c) 2024 CERN for the benefit of the FastCaloSim project

#ifndef TFCSHistoLateralShapeWeight_h
#define TFCSHistoLateralShapeWeight_h

#include "FastCaloSim/Core/TFCSLateralShapeParametrizationHitBase.h"

class TH1;
namespace FastCaloSim::Core
{
class TFCSHistoLateralShapeWeight
    : public TFCSLateralShapeParametrizationHitBase
{
public:
  TFCSHistoLateralShapeWeight(const char* name = nullptr,
                              const char* title = nullptr);
  virtual ~TFCSHistoLateralShapeWeight();

  /// weight the energy of one hit in order to generate fluctuations
  virtual FCSReturnCode simulate_hit(
      Hit& hit,
      TFCSSimulationState& simulstate,
      const TFCSTruthState* truth,
      const TFCSExtrapolationState* extrapol) override;

  /// Init from histogram. The integral of the histogram is used as number of
  /// expected hits to be generated
  bool Initialize(TH1* hist);

  TH1* getHistogram() const { return m_hist; };

  virtual void Print(Option_t* option = "") const override;
  virtual void setMinWeight(float minWeight) { m_minWeight = minWeight; }
  virtual void setMaxWeight(float maxWeight) { m_maxWeight = maxWeight; }
  virtual float getMinWeight() const override;
  virtual float getMaxWeight() const override;

protected:
  /// Histogram to be used for the shape simulation
  /// The histogram x-axis should be in dR^2=deta^2+dphi^2
  TH1* m_hist {nullptr};
  float m_minWeight {-1.};
  float m_maxWeight {-1.};

  ClassDefOverride(TFCSHistoLateralShapeWeight,
                   2)  // TFCSHistoLateralShapeWeight
};
}  // namespace FastCaloSim::Core
#endif
