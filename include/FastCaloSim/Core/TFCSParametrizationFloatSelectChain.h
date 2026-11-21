// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#pragma once

#include "FastCaloSim/Core/TFCSParametrizationBinnedChain.h"

// Custom copy constructor, moves fall back to it
// NOLINTNEXTLINE(cppcoreguidelines-special-member-functions,hicpp-special-member-functions)
class TFCSParametrizationFloatSelectChain
    : public TFCSParametrizationBinnedChain
{
public:
  explicit TFCSParametrizationFloatSelectChain(const char* name = nullptr,
                                               const char* title = nullptr)
      : TFCSParametrizationBinnedChain(name, title)
      , m_bin_low_edge(1, 0)
  {
  }
  TFCSParametrizationFloatSelectChain(
      const TFCSParametrizationFloatSelectChain& ref)
      : TFCSParametrizationBinnedChain(ref)
      , m_bin_low_edge(1, 0)
  {
  }

  virtual auto push_back_in_bin(TFCSParametrizationBase* param,
                                float low,
                                float up) -> std::optional<std::size_t>;
  /// Should not be used unless the bin boundaries are already defined!
  void push_back_in_bin(TFCSParametrizationBase* param,
                        std::size_t bin) override;

  // returns no value if outside range
  auto val_to_bin(float val) const -> std::optional<std::size_t>;

  virtual auto get_bin_low_edge(std::size_t bin) const -> double
  {
    return static_cast<double>(m_bin_low_edge[bin]);
  }
  virtual auto get_bin_up_edge(std::size_t bin) const -> double
  {
    return static_cast<double>(m_bin_low_edge[bin + 1]);
  }

protected:
  std::vector<float> m_bin_low_edge;

private:
  ClassDefOverride(TFCSParametrizationFloatSelectChain,
                   1)  // TFCSParametrizationFloatSelectChain
};
