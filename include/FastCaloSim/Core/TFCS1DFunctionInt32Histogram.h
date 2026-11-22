// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#pragma once

#include <cstdint>
#include <utility>
#include <vector>

#include "FastCaloSim/Core/TFCS1DFunction.h"
#include "FastCaloSim/Definitions/GuideTable.h"

class TH2;

class TFCS1DFunctionInt32Histogram : public TFCS1DFunction
{
public:
  explicit TFCS1DFunctionInt32Histogram(const TH1* hist = nullptr)
  {
    if (hist != nullptr) {
      Initialize(hist);
    }
  }
  ~TFCS1DFunctionInt32Histogram() override = default;
  TFCS1DFunctionInt32Histogram(const TFCS1DFunctionInt32Histogram&) = default;
  TFCS1DFunctionInt32Histogram(TFCS1DFunctionInt32Histogram&&) = default;
  auto operator=(const TFCS1DFunctionInt32Histogram&)
      -> TFCS1DFunctionInt32Histogram& = default;
  auto operator=(TFCS1DFunctionInt32Histogram&&)
      -> TFCS1DFunctionInt32Histogram& = default;

  void Initialize(const TH1* hist);

  using TFCS1DFunction::rnd_to_fct;

  using HistoContent_t = std::uint32_t;
  static const HistoContent_t s_MaxValue;

  /// Function gets random number rnd in the range [0,1) as argument
  /// and returns function value according to a histogram distribution
  auto rnd_to_fct(double rnd) const -> double override;

  auto operator==(const TFCS1DFunction& ref) const -> bool override;

  auto get_HistoBordersx() const -> const std::vector<float>&
  {
    return m_HistoBorders;
  }
  auto get_HistoBordersx() -> std::vector<float>& { return m_HistoBorders; }
  auto get_HistoContents() const -> const std::vector<HistoContent_t>&
  {
    return m_HistoContents;
  }
  /// Replaces the cumulative contents (and rebuilds the search table)
  void set_HistoContents(std::vector<HistoContent_t> contents)
  {
    m_HistoContents = std::move(contents);
    m_guide.build(m_HistoContents);
  }
  /// Search table for m_HistoContents, built from it after reading
  auto guide_table() const -> const FastCaloSim::GuideTable<HistoContent_t>&
  {
    return m_guide;
  }

protected:
  std::vector<float> m_HistoBorders;
  std::vector<HistoContent_t> m_HistoContents;
  /// Transient, built from m_HistoContents (also by a read rule after reading)
  FastCaloSim::GuideTable<HistoContent_t> m_guide;  //!

private:
  ClassDefOverride(TFCS1DFunctionInt32Histogram,
                   1)  // TFCS1DFunctionInt32Histogram
};
