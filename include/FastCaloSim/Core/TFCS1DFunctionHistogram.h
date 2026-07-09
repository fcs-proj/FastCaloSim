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
  TFCS1DFunctionHistogram(TH1* hist, double);

  void Initialize(TH1* hist, double);

  using TFCS1DFunction::rnd_to_fct;
  auto rnd_to_fct(double rnd) const -> double override;
  auto vector_to_histo() -> TH1*;
  auto get_inverse(double rnd) const -> double;
  static auto linear(double y1, double y2, double x1, double x2, double x)
      -> double;
  static auto non_linear(double y1, double y2, double x1, double x2, double x)
      -> double;

  static auto get_maxdev(TH1*, TH1D*) -> double;
  void smart_rebin_loop(TH1* hist, double);
  static auto get_change(TH1*) -> double;
  static auto smart_rebin(TH1D*) -> TH1D*;
  static auto histo_to_array(TH1*) -> std::unique_ptr<double[]>;
  static auto sample_from_histo(TH1* hist, double) -> double;
  auto sample_from_histovalues(double) -> double;

  auto get_HistoBorders() const -> const std::vector<float>&
  {
    return m_HistoBorders;
  };
  auto get_HistoContents() const -> const std::vector<float>&
  {
    return m_HistoContents;
  };

protected:
  std::vector<float> m_HistoBorders;
  std::vector<float> m_HistoContents;

  ClassDefOverride(TFCS1DFunctionHistogram, 1)  // TFCS1DFunctionHistogram
};
