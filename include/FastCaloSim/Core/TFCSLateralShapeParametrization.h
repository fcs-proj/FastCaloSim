// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#pragma once

#include <cstddef>
#include <optional>
#include <utility>

#include "FastCaloSim/Core/TFCSParametrization.h"

class TFCSLateralShapeParametrization : public TFCSParametrization
{
public:
  explicit TFCSLateralShapeParametrization(const char* name = nullptr,
                                           const char* title = nullptr);

  auto is_match_Ekin_bin(std::size_t bin) const -> bool override
  {
    return !Ekin_bin() || bin == *Ekin_bin();
  }
  auto is_match_calosample(std::size_t calosample) const -> bool override
  {
    return std::cmp_equal(m_calosample, calosample);
  }

  auto is_match_all_Ekin_bin() const -> bool override { return !Ekin_bin(); }
  auto is_match_all_calosample() const -> bool override { return false; }

  /// The Ekin bin; no value for all bins
  auto Ekin_bin() const -> std::optional<std::size_t>
  {
    if (m_Ekin_bin < 0) {
      return std::nullopt;
    }
    return static_cast<std::size_t>(m_Ekin_bin);
  }
  void set_Ekin_bin(std::size_t bin);

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
  // Not virtual: compares only the members of this class, the derived classes
  // call it qualified
  // NOLINTNEXTLINE(bugprone-derived-method-shadowing-base-method): by design
  auto compare(const TFCSParametrizationBase& ref) const -> bool;

private:
  // Persistent as int, with -1 for all bins
  int m_Ekin_bin {-1};
  // Persistent as int, with -1 for no layer
  int m_calosample {-1};

  ClassDefOverride(TFCSLateralShapeParametrization,
                   1)  // TFCSLateralShapeParametrization
};
