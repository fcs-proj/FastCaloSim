// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#pragma once

#include <vector>

#include "FastCaloSim/Core/TFCS1DFunction.h"

class TH2;

class TFCS1DFunctionInt32Histogram : public TFCS1DFunction
{
public:
  TFCS1DFunctionInt32Histogram(const TH1* hist = nullptr)
  {
    if (hist)
      Initialize(hist);
  };
  ~TFCS1DFunctionInt32Histogram() override = default;

  void Initialize(const TH1* hist);

  using TFCS1DFunction::rnd_to_fct;

  using HistoContent_t = uint32_t;
  static const HistoContent_t s_MaxValue;

  /// Function gets random number rnd in the range [0,1) as argument
  /// and returns function value according to a histogram distribution
  auto rnd_to_fct(double rnd) const -> double override;

  auto operator==(const TFCS1DFunction& ref) const -> bool override;

  auto get_HistoBordersx() const -> const std::vector<float>&
  {
    return m_HistoBorders;
  };
  auto get_HistoBordersx() -> std::vector<float>& { return m_HistoBorders; };
  auto get_HistoContents() const -> const std::vector<HistoContent_t>&
  {
    return m_HistoContents;
  };
  auto get_HistoContents() -> std::vector<HistoContent_t>&
  {
    return m_HistoContents;
  };

protected:
  std::vector<float> m_HistoBorders;
  std::vector<HistoContent_t> m_HistoContents;

private:
  ClassDef(TFCS1DFunctionInt32Histogram, 1)  // TFCS1DFunctionInt32Histogram
};
