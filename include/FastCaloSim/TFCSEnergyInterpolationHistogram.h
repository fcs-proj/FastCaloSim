/*
  Copyright (C) 2002-2020 CERN for the benefit of the ATLAS collaboration
*/

#ifndef ISF_FASTCALOSIMEVENT_TFCSEnergyInterpolationHistogram_h
#define ISF_FASTCALOSIMEVENT_TFCSEnergyInterpolationHistogram_h

#include "FastCaloSim/TFCSParametrization.h"
#include "TH1F.h"

class TFCSEnergyInterpolationHistogram : public TFCSParametrization {
public:
  TFCSEnergyInterpolationHistogram(const char *name = nullptr,
                                   const char *title = nullptr);

  enum FCSEnergyInitializationStatusBits {
    kOnlyScaleEnergy =
        BIT(15) ///< Set this bit in the TObject bit field the simulated energy
                ///< should only be scaled by the spline
  };

  bool OnlyScaleEnergy() const { return TestBit(kOnlyScaleEnergy); };
  void set_OnlyScaleEnergy() { SetBit(kOnlyScaleEnergy); };
  void reset_OnlyScaleEnergy() { ResetBit(kOnlyScaleEnergy); };

  virtual bool is_match_Ekin_bin(int /*Ekin_bin*/) const override {
    return true;
  };
  virtual bool is_match_calosample(int /*calosample*/) const override {
    return true;
  };

  /// Initialize interpolation from histogram
  /// x values should be Ekin, y values should <E(reco)/Ekin(true)>
  void InitFromHist(const TH1F &hist) { m_hist = hist; };

  const TH1F &hist() const { return m_hist; };

  /// Initialize simulstate with the mean reconstructed energy in the
  /// calorimater expeted from the true kinetic energy
  virtual FCSReturnCode
  simulate(TFCSSimulationState &simulstate, const TFCSTruthState *truth,
           const TFCSExtrapolationState *extrapol) const override;
  void Print(Option_t *option = "") const override;

  static void unit_test(TFCSSimulationState *simulstate = nullptr,
                        TFCSTruthState *truth = nullptr,
                        const TFCSExtrapolationState *extrapol = nullptr,
                        TH1F *hist = nullptr);

private:
  TH1F m_hist;

  ClassDefOverride(TFCSEnergyInterpolationHistogram,
                   1) // TFCSEnergyInterpolationHistogram
};

#if defined(__ROOTCLING__) && defined(__FastCaloSimStandAlone__)
#pragma link C++ class TFCSEnergyInterpolationHistogram + ;
#endif

#endif
