// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#ifndef ISF_FASTCALOSIMEVENT_TFCSEnergyInterpolationHistogram_h
#define ISF_FASTCALOSIMEVENT_TFCSEnergyInterpolationHistogram_h

#include "FastCaloSim/Core/TFCSParametrization.h"
#include "TH1F.h"

class TFCSEnergyInterpolationHistogram : public TFCSParametrization
{
public:
  TFCSEnergyInterpolationHistogram(const char* name = nullptr,
                                   const char* title = nullptr);

  enum FCSEnergyInitializationStatusBits
  {
    kOnlyScaleEnergy =
        BIT(15)  ///< Set this bit in the TObject bit field the simulated energy
                 ///< should only be scaled by the spline
  };

  auto OnlyScaleEnergy() const -> bool { return TestBit(kOnlyScaleEnergy); };
  void set_OnlyScaleEnergy() { SetBit(kOnlyScaleEnergy); };
  void reset_OnlyScaleEnergy() { ResetBit(kOnlyScaleEnergy); };

  auto is_match_Ekin_bin(int /*Ekin_bin*/) const -> bool override
  {
    return true;
  };
  auto is_match_calosample(int /*calosample*/) const -> bool override
  {
    return true;
  };

  /// Initialize interpolation from histogram
  /// x values should be Ekin, y values should <E(reco)/Ekin(true)>
  void InitFromHist(const TH1F& hist) { m_hist = hist; };

  auto hist() const -> const TH1F& { return m_hist; };

  /// Initialize simulstate with the mean reconstructed energy in the
  /// calorimater expected from the true kinetic energy
  auto simulate(TFCSSimulationState& simulstate,
                const TFCSTruthState* truth,
                const TFCSExtrapolationState* extrapol) const
      -> FCSReturnCode override;
  void Print(Option_t* option = "") const override;

private:
  TH1F m_hist;

  ClassDefOverride(TFCSEnergyInterpolationHistogram,
                   1)  // TFCSEnergyInterpolationHistogram
};

#if defined(__ROOTCLING__) && defined(__FastCaloSimStandAlone__)
#  pragma link C++ class TFCSEnergyInterpolationHistogram + ;
#endif

#endif
