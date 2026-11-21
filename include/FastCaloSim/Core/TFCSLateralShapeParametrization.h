// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#pragma once

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
  }
  auto is_match_calosample(std::size_t calosample) const -> bool override
  {
    return m_calosample >= 0
        && static_cast<std::size_t>(m_calosample) == calosample;
  }

  auto is_match_all_Ekin_bin() const -> bool override
  {
    if (Ekin_bin() == -1)
      return true;
    return false;
  }
  auto is_match_all_calosample() const -> bool override { return false; }

  auto Ekin_bin() const -> int { return m_Ekin_bin; }
  void set_Ekin_bin(int bin);

  /// The calorimeter layer; only valid if has_calosample()
  auto calosample() const -> std::size_t
  {
    return static_cast<std::size_t>(m_calosample);
  }
  auto has_calosample() const -> bool { return m_calosample >= 0; }
  void set_calosample(std::size_t cs);

  virtual void set_pdgid_Ekin_eta_Ekin_bin_calosample(
      const TFCSLateralShapeParametrization& ref);

  void Print(Option_t* option = "") const override;

protected:
  auto compare(const TFCSParametrizationBase& ref) const -> bool;

private:
  int m_Ekin_bin;
  // Persistent as int, with -1 for no layer
  int m_calosample;

  ClassDefOverride(TFCSLateralShapeParametrization,
                   1)  // TFCSLateralShapeParametrization
};
