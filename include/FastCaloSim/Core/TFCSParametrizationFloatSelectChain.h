// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#pragma once

#include "FastCaloSim/Core/TFCSParametrizationBinnedChain.h"

class TFCSParametrizationFloatSelectChain
    : public TFCSParametrizationBinnedChain
{
public:
  TFCSParametrizationFloatSelectChain(const char* name = nullptr,
                                      const char* title = nullptr)
      : TFCSParametrizationBinnedChain(name, title)
      , m_bin_low_edge(1, 0) {};
  TFCSParametrizationFloatSelectChain(
      const TFCSParametrizationFloatSelectChain& ref)
      : TFCSParametrizationBinnedChain(ref)
      , m_bin_low_edge(1, 0) {};

  virtual auto push_back_in_bin(TFCSParametrizationBase* param,
                                float low,
                                float up) -> int;
  /// Should not be used unless the bin boundaries are already defined!
  void push_back_in_bin(TFCSParametrizationBase* param,
                        unsigned int bin) override;

  // return -1 if outside range
  auto val_to_bin(float val) const -> int;

  virtual auto get_bin_low_edge(int bin) const -> double
  {
    return m_bin_low_edge[bin];
  };
  virtual auto get_bin_up_edge(int bin) const -> double
  {
    return m_bin_low_edge[bin + 1];
  };

protected:
  std::vector<float> m_bin_low_edge;

private:
  ClassDefOverride(TFCSParametrizationFloatSelectChain,
                   1)  // TFCSParametrizationFloatSelectChain
};
