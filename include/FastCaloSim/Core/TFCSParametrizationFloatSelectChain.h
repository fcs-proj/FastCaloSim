// Copyright (c) 2024 CERN for the benefit of the FastCaloSim project

#ifndef TFCSParametrizationFloatSelectChain_h
#define TFCSParametrizationFloatSelectChain_h

#include "FastCaloSim/Core/TFCSParametrizationBinnedChain.h"

namespace FastCaloSim::Core
{

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

  virtual int push_back_in_bin(TFCSParametrizationBase* param,
                               float low,
                               float up);
  /// Should not be used unless the bin boundaries are already defined!
  virtual void push_back_in_bin(TFCSParametrizationBase* param,
                                unsigned int bin) override;

  // return -1 if outside range
  int val_to_bin(float val) const;

  virtual double get_bin_low_edge(int bin) const
  {
    return m_bin_low_edge[bin];
  };
  virtual double get_bin_up_edge(int bin) const
  {
    return m_bin_low_edge[bin + 1];
  };

protected:
  std::vector<float> m_bin_low_edge;

private:
  ClassDefOverride(TFCSParametrizationFloatSelectChain,
                   1)  // TFCSParametrizationFloatSelectChain
};
}  // namespace FastCaloSim::Core
#endif
