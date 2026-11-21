// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#pragma once

#include <memory>
#include <vector>

#include "FastCaloSim/Core/TFCS1DFunction.h"
#include "TH1.h"

class TFCS1DFunctionHistogram : public TFCS1DFunction
{
public:
  TFCS1DFunctionHistogram() = default;
  TFCS1DFunctionHistogram(TH1* hist, double cut_maxdev);

  void Initialize(TH1* hist, double cut_maxdev);

  using TFCS1DFunction::rnd_to_fct;
  auto rnd_to_fct(double rnd) const -> double override;
  auto get_inverse(double rnd) const -> double;
  static auto non_linear(double y1, double y2, double x1, double x2, double x)
      -> double;

  static auto get_maxdev(TH1*, TH1D*) -> double;
  void smart_rebin_loop(TH1* hist, double);
  static auto get_change(TH1*) -> double;
  static auto smart_rebin(TH1D*) -> std::unique_ptr<TH1D>;

  auto get_HistoBorders() const -> const std::vector<float>&
  {
    return m_HistoBorders;
  }
  auto get_HistoContents() const -> const std::vector<float>&
  {
    return m_HistoContents;
  }

protected:
  std::vector<float> m_HistoBorders;
  std::vector<float> m_HistoContents;

  ClassDefOverride(TFCS1DFunctionHistogram, 1)  // TFCS1DFunctionHistogram
};
