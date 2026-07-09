// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#ifndef TFCSHistoLateralShapeWeight_h
#define TFCSHistoLateralShapeWeight_h

#include "FastCaloSim/Core/TFCSLateralShapeParametrizationHitBase.h"

class TH1;

class TFCSHistoLateralShapeWeight
    : public TFCSLateralShapeParametrizationHitBase
{
public:
  TFCSHistoLateralShapeWeight(const char* name = nullptr,
                              const char* title = nullptr);
  ~TFCSHistoLateralShapeWeight() override;

  /// weight the energy of one hit in order to generate fluctuations
  auto simulate_hit(Hit& hit,
                    TFCSSimulationState& simulstate,
                    const TFCSTruthState* truth,
                    const TFCSExtrapolationState* extrapol)
      -> FCSReturnCode override;

  /// Init from histogram. The integral of the histogram is used as number of
  /// expected hits to be generated
  auto Initialize(TH1* hist) -> bool;

  auto getHistogram() const -> TH1* { return m_hist; };

  void Print(Option_t* option = "") const override;
  virtual void setMinWeight(float minWeight) { m_minWeight = minWeight; }
  virtual void setMaxWeight(float maxWeight) { m_maxWeight = maxWeight; }
  auto getMinWeight() const -> float override;
  auto getMaxWeight() const -> float override;

protected:
  /// Histogram to be used for the shape simulation
  /// The histogram x-axis should be in dR^2=deta^2+dphi^2
  TH1* m_hist {nullptr};
  float m_minWeight {-1.};
  float m_maxWeight {-1.};

  ClassDefOverride(TFCSHistoLateralShapeWeight,
                   2)  // TFCSHistoLateralShapeWeight
};

#endif
