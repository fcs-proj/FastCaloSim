// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#pragma once

#include <vector>

#include "FastCaloSim/Core/TFCS2DFunction.h"

class TH2;

class TFCS2DFunctionHistogram : public TFCS2DFunction
{
public:
  explicit TFCS2DFunctionHistogram(TH2* hist = nullptr)
  {
    if (hist != nullptr) {
      Initialize(hist);
    }
  }
  ~TFCS2DFunctionHistogram() override = default;
  TFCS2DFunctionHistogram(const TFCS2DFunctionHistogram&) = default;
  TFCS2DFunctionHistogram(TFCS2DFunctionHistogram&&) = default;
  auto operator=(const TFCS2DFunctionHistogram&)
      -> TFCS2DFunctionHistogram& = default;
  auto operator=(TFCS2DFunctionHistogram&&)
      -> TFCS2DFunctionHistogram& = default;

  void Initialize(TH2* hist);

  using TFCS2DFunction::rnd_to_fct;
  void rnd_to_fct(double& valuex,
                  double& valuey,
                  double rnd0,
                  double rnd1) const override;

  auto get_HistoBordersx() const -> const std::vector<float>&
  {
    return m_HistoBorders;
  }
  auto get_HistoBordersx() -> std::vector<float>& { return m_HistoBorders; }
  auto get_HistoBordersy() const -> const std::vector<float>&
  {
    return m_HistoBordersy;
  }
  auto get_HistoBordersy() -> std::vector<float>& { return m_HistoBordersy; }
  auto get_HistoContents() const -> const std::vector<float>&
  {
    return m_HistoContents;
  }
  auto get_HistoContents() -> std::vector<float>& { return m_HistoContents; }

protected:
  std::vector<float> m_HistoBorders;
  std::vector<float> m_HistoBordersy;
  std::vector<float> m_HistoContents;

private:
  ClassDefOverride(TFCS2DFunctionHistogram, 1)  // TFCS2DFunctionHistogram
};
