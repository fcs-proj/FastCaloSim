// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#pragma once

#include "FastCaloSim/Core/TFCSParametrizationBinnedChain.h"

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
      const TFCSParametrizationFloatSelectChain& ref) = default;
  TFCSParametrizationFloatSelectChain(
      TFCSParametrizationFloatSelectChain&& ref) = default;
  auto operator=(const TFCSParametrizationFloatSelectChain& ref)
      -> TFCSParametrizationFloatSelectChain& = default;
  auto operator=(TFCSParametrizationFloatSelectChain&& ref)
      -> TFCSParametrizationFloatSelectChain& = default;
  ~TFCSParametrizationFloatSelectChain() override = default;

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
