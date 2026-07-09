// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#ifndef TFCSLateralShapeParametrization_h
#define TFCSLateralShapeParametrization_h

#include "FastCaloSim/Core/TFCSParametrization.h"

class TFCSLateralShapeParametrization : public TFCSParametrization
{
public:
  TFCSLateralShapeParametrization(const char* name = nullptr,
                                  const char* title = nullptr);

  auto is_match_Ekin_bin(int bin) const -> bool override
  {
    if (Ekin_bin() == -1)
      return true;
    return bin == Ekin_bin();
  };
  auto is_match_calosample(int calosample) const -> bool override
  {
    return calosample == m_calosample;
  };

  auto is_match_all_Ekin_bin() const -> bool override
  {
    if (Ekin_bin() == -1)
      return true;
    return false;
  };
  auto is_match_all_calosample() const -> bool override { return false; };

  auto Ekin_bin() const -> int { return m_Ekin_bin; };
  void set_Ekin_bin(int bin);

  auto calosample() const -> int { return m_calosample; };
  void set_calosample(int cs);

  virtual void set_pdgid_Ekin_eta_Ekin_bin_calosample(
      const TFCSLateralShapeParametrization& ref);

  void Print(Option_t* option = "") const override;

protected:
  auto compare(const TFCSParametrizationBase& ref) const -> bool;

private:
  int m_Ekin_bin;
  int m_calosample;

  ClassDefOverride(TFCSLateralShapeParametrization,
                   1)  // TFCSLateralShapeParametrization
};

#endif
