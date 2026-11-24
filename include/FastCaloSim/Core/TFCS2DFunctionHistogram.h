// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#pragma once

#include <utility>
#include <vector>

#include "FastCaloSim/Core/TFCS2DFunction.h"
#include "FastCaloSim/Definitions/GuideTable.h"

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

  /// Const or not as the object
  auto get_HistoBordersx(this auto&& self) -> auto&&
  {
    return std::forward_like<decltype(self)>(self.m_HistoBorders);
  }
  /// Const or not as the object
  auto get_HistoBordersy(this auto&& self) -> auto&&
  {
    return std::forward_like<decltype(self)>(self.m_HistoBordersy);
  }
  auto get_HistoContents() const -> const std::vector<float>&
  {
    return m_HistoContents;
  }
  /// Replaces the cumulative contents (and rebuilds the search table)
  void set_HistoContents(std::vector<float> contents)
  {
    m_HistoContents = std::move(contents);
    m_guide.build(m_HistoContents);
  }
  /// Search table for m_HistoContents, built from it after reading
  auto guide_table() const -> const FastCaloSim::GuideTable<float>&
  {
    return m_guide;
  }

protected:
  std::vector<float> m_HistoBorders;
  std::vector<float> m_HistoBordersy;
  std::vector<float> m_HistoContents;
  /// Transient, built from m_HistoContents (also by a read rule after reading)
  FastCaloSim::GuideTable<float> m_guide;  //!

private:
  ClassDefOverride(TFCS2DFunctionHistogram, 1)  // TFCS2DFunctionHistogram
};
