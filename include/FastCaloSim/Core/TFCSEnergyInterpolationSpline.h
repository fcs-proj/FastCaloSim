// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#ifndef ISF_FASTCALOSIMEVENT_TFCSEnergyInterpolationSpline_h
#define ISF_FASTCALOSIMEVENT_TFCSEnergyInterpolationSpline_h

#include "FastCaloSim/Core/TFCSParametrization.h"
#include "TSpline.h"

class TGraph;

class TFCSEnergyInterpolationSpline : public TFCSParametrization
{
public:
  TFCSEnergyInterpolationSpline(const char* name = nullptr,
                                const char* title = nullptr);

  /// Status bit for energy initialization
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

  /// Initialize interpolation from spline
  /// x values should be log(Ekin), y values should <E(reco)/Ekin(true)>
  void InitFromSpline(const TSpline3& spline) { m_spline = spline; };

  /// Initialize spline interpolation from arrays in log(Ekin) and
  /// response=<E(reco)/Ekin(true)> opt, valbeg and valend as defined for
  /// TSpline3
  void InitFromArrayInLogEkin(Int_t np,
                              Double_t logEkin[],
                              Double_t response[],
                              const char* opt = nullptr,
                              Double_t valbeg = 0,
                              Double_t valend = 0);

  /// Initialize spline interpolation from arrays in Ekin and
  /// response=<E(reco)/Ekin(true)> opt, valbeg and valend as defined for
  /// TSpline3
  void InitFromArrayInEkin(Int_t np,
                           Double_t Ekin[],
                           Double_t response[],
                           const char* opt = nullptr,
                           Double_t valbeg = 0,
                           Double_t valend = 0);

  auto spline() const -> const TSpline3& { return m_spline; };

  /// Initialize simulstate with the mean reconstructed energy in the
  /// calorimater expected from the true kinetic energy
  auto simulate(TFCSSimulationState& simulstate,
                const TFCSTruthState* truth,
                const TFCSExtrapolationState* extrapol) const
      -> FCSReturnCode override;

  void Print(Option_t* option = "") const override;

private:
  TSpline3 m_spline;

  ClassDefOverride(TFCSEnergyInterpolationSpline,
                   1)  // TFCSEnergyInterpolationSpline
};

#endif
