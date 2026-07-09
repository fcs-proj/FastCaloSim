// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#ifndef TFCSHistoLateralShapeParametrization_h
#define TFCSHistoLateralShapeParametrization_h

#include "FastCaloSim/Core/TFCS2DFunctionHistogram.h"
#include "FastCaloSim/Core/TFCSLateralShapeParametrizationHitBase.h"
#include "FastCaloSim/Core/TFCSTruthState.h"

class TH2;

class TFCSHistoLateralShapeParametrization
    : public TFCSLateralShapeParametrizationHitBase
{
public:
  TFCSHistoLateralShapeParametrization(const char* name = nullptr,
                                       const char* title = nullptr);
  ~TFCSHistoLateralShapeParametrization() override;

  /// Status bit for FCS needs
  enum FCSStatusBits
  {
    k_phi_symmetric =
        BIT(15)  ///< Set this bit to simulate phi symmetric histograms
  };

  /// will actually not store the geometry information, but rather used to check
  /// the validity of the 2D shape histogram
  void set_geometry(CaloGeo* geo) override;

  auto is_phi_symmetric() const -> bool { return TestBit(k_phi_symmetric); };
  virtual void set_phi_symmetric() { SetBit(k_phi_symmetric); };
  virtual void reset_phi_symmetric() { ResetBit(k_phi_symmetric); };

  /// set the integral of the histogram to the desired number of hits
  void set_number_of_hits(float nhits);

  auto get_number_of_expected_hits() const -> float { return m_nhits; };

  /// set an offset in r on the simulated histogram
  void set_r_offset(float r_offset) { m_r_offset = r_offset; };
  auto r_offset() const -> float { return m_r_offset; };

  /// set an scale factor for r on the simulated histogram
  void set_r_scale(float r_scale) { m_r_scale = r_scale; };
  auto r_scale() const -> float { return m_r_scale; };

  /// default for this class is to simulate get_number_of_expected_hits() hits,
  /// which gives fluctuations sigma^2=1/get_number_of_expected_hits()
  auto get_sigma2_fluctuation(TFCSSimulationState& simulstate,
                              const TFCSTruthState* truth,
                              const TFCSExtrapolationState* extrapol) const
      -> double override;

  /// default for this class is to simulate get_number_of_expected_hits() hits
  auto get_number_of_hits(TFCSSimulationState& simulstate,
                          const TFCSTruthState* truth,
                          const TFCSExtrapolationState* extrapol) const
      -> int override;

  /// simulated one hit position with weight that should be put into simulstate
  /// sometime later all hit weights should be resacled such that their final
  /// sum is simulstate->E(sample) someone also needs to map all hits into cells
  auto simulate_hit(Hit& hit,
                    TFCSSimulationState& simulstate,
                    const TFCSTruthState* truth,
                    const TFCSExtrapolationState* extrapol)
      -> FCSReturnCode override;

  /// Init from histogram. The integral of the histogram is used as number of
  /// expected hits to be generated
  auto Initialize(TH2* hist) -> bool;
  auto Initialize(const char* filepath, const char* histname) -> bool;

  auto histogram() -> TFCS2DFunctionHistogram& { return m_hist; };
  auto histogram() const -> const TFCS2DFunctionHistogram& { return m_hist; };

  void Print(Option_t* option = "") const override;

protected:
  /// Histogram to be used for the shape simulation
  TFCS2DFunctionHistogram m_hist;
  float m_nhits;
  float m_r_offset;
  float m_r_scale;

private:
  ClassDefOverride(TFCSHistoLateralShapeParametrization,
                   2)  // TFCSHistoLateralShapeParametrization
};

#endif
