// Copyright (c) 2024 CERN for the benefit of the FastCaloSim project

#ifndef TFCS1DFunctionInt32Histogram_h
#define TFCS1DFunctionInt32Histogram_h

#include <vector>

#include "FastCaloSim/Core/TFCS1DFunction.h"

class TH2;

namespace FastCaloSim::Core
{

class TFCS1DFunctionInt32Histogram : public TFCS1DFunction
{
public:
  TFCS1DFunctionInt32Histogram(const TH1* hist = nullptr)
  {
    if (hist)
      Initialize(hist);
  };
  ~TFCS1DFunctionInt32Histogram() {};

  void Initialize(const TH1* hist);

  using TFCS1DFunction::rnd_to_fct;

  typedef uint32_t HistoContent_t;
  static const HistoContent_t s_MaxValue;

  /// Function gets random number rnd in the range [0,1) as argument
  /// and returns function value according to a histogram distribution
  virtual double rnd_to_fct(double rnd) const;

  virtual bool operator==(const TFCS1DFunction& ref) const;

  const std::vector<float>& get_HistoBordersx() const
  {
    return m_HistoBorders;
  };
  std::vector<float>& get_HistoBordersx() { return m_HistoBorders; };
  const std::vector<HistoContent_t>& get_HistoContents() const
  {
    return m_HistoContents;
  };
  std::vector<HistoContent_t>& get_HistoContents() { return m_HistoContents; };

protected:
  std::vector<float> m_HistoBorders;
  std::vector<HistoContent_t> m_HistoContents;

private:
  ClassDef(TFCS1DFunctionInt32Histogram, 1)  // TFCS1DFunctionInt32Histogram
};
}  // namespace FastCaloSim::Core

#endif
