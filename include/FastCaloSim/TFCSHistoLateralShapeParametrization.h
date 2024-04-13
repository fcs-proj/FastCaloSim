/*
  Copyright (C) 2002-2019 CERN for the benefit of the ATLAS collaboration
*/

#ifndef TFCSHistoLateralShapeParametrization_h
#define TFCSHistoLateralShapeParametrization_h

#include "FastCaloSim/TFCSLateralShapeParametrizationHitBase.h"
#include "FastCaloSim/TFCS2DFunctionHistogram.h"
#include "FastCaloSim/TFCSTruthState.h"
#ifdef USE_GPU
#include "ISF_FastCaloGpu/LoadGpuFuncHist.h"
#include "ISF_FastCaloGpu/FH_structs.h"
#endif

class TH2;
class ICaloGeometry;

class TFCSHistoLateralShapeParametrization
    : public TFCSLateralShapeParametrizationHitBase {
public:
  TFCSHistoLateralShapeParametrization(const char *name = nullptr,
                                       const char *title = nullptr);
  ~TFCSHistoLateralShapeParametrization();

  /// Status bit for FCS needs
  enum FCSStatusBits {
    k_phi_symmetric =
        BIT(15) ///< Set this bit to simulate phi symmetric histograms
  };

  /// will actually not store the geometry information, but rather used to check
  /// the validity of the 2D shape histogram
  virtual void set_geometry(ICaloGeometry *geo) override;

  bool is_phi_symmetric() const { return TestBit(k_phi_symmetric); };
  virtual void set_phi_symmetric() { SetBit(k_phi_symmetric); };
  virtual void reset_phi_symmetric() { ResetBit(k_phi_symmetric); };

  /// set the integral of the histogram to the desired number of hits
  void set_number_of_hits(float nhits);

  float get_number_of_expected_hits() const { return m_nhits; };

  /// set an offset in r on the simulated histogram
  void set_r_offset(float r_offset) { m_r_offset = r_offset; };
  float r_offset() const { return m_r_offset; };

  /// set an scale factor for r on the simulated histogram
  void set_r_scale(float r_scale) { m_r_scale = r_scale; };
  float r_scale() const { return m_r_scale; };

  /// default for this class is to simulate get_number_of_expected_hits() hits,
  /// which gives fluctuations sigma^2=1/get_number_of_expected_hits()
  virtual double
  get_sigma2_fluctuation(TFCSSimulationState &simulstate,
                         const TFCSTruthState *truth,
                         const TFCSExtrapolationState *extrapol) const override;

  /// default for this class is to simulate get_number_of_expected_hits() hits
  virtual int
  get_number_of_hits(TFCSSimulationState &simulstate,
                     const TFCSTruthState *truth,
                     const TFCSExtrapolationState *extrapol) const override;

  /// simulated one hit position with weight that should be put into simulstate
  /// sometime later all hit weights should be resacled such that their final
  /// sum is simulstate->E(sample) someone also needs to map all hits into cells
  virtual FCSReturnCode
  simulate_hit(Hit &hit, TFCSSimulationState &simulstate,
               const TFCSTruthState *truth,
               const TFCSExtrapolationState *extrapol) override;

  /// Init from histogram. The integral of the histogram is used as number of
  /// expected hits to be generated
  bool Initialize(TH2 *hist);
  bool Initialize(const char *filepath, const char *histname);

  TFCS2DFunctionHistogram &histogram() { return m_hist; };
  const TFCS2DFunctionHistogram &histogram() const { return m_hist; };

  void Print(Option_t *option = "") const override;

#ifdef USE_GPU
  // for FCS-GPU
  // construct the hist function and copy to GPU
  // will not compile by default
  void set_d_HistFunc(FH2D *hf_ptr) { m_d_HistFunc = hf_ptr; };
  const FH2D *d_HistFunc() { return m_d_HistFunc; };
  LoadGpuFuncHist *LdFH() { return m_LdFH; };
  void LoadHistFuncs();
#endif

protected:
  /// Histogram to be used for the shape simulation
  TFCS2DFunctionHistogram m_hist;
  float m_nhits;
  float m_r_offset;
  float m_r_scale;

private:
#ifdef USE_GPU
  FH2D *m_d_HistFunc = nullptr;
  LoadGpuFuncHist *m_LdFH = nullptr;
#endif

  ClassDefOverride(TFCSHistoLateralShapeParametrization,
                   2) // TFCSHistoLateralShapeParametrization
};

#endif
