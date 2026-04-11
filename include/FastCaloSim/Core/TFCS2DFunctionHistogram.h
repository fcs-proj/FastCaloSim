// Copyright (c) 2024 CERN for the benefit of the FastCaloSim project

#ifndef TFCS2DFunctionHistogram_h
#define TFCS2DFunctionHistogram_h

#include <vector>

#include "FastCaloSim/Core/TFCS2DFunction.h"

class TH2;

namespace FastCaloSim::Core
{
class TFCS2DFunctionHistogram : public TFCS2DFunction
{
public:
  TFCS2DFunctionHistogram(TH2* hist = nullptr)
  {
    if (hist)
      Initialize(hist);
  };
  ~TFCS2DFunctionHistogram() {};

  void Initialize(TH2* hist);

  using TFCS2DFunction::rnd_to_fct;
  virtual void rnd_to_fct(float& valuex,
                          float& valuey,
                          float rnd0,
                          float rnd1) const;

  const std::vector<float>& get_HistoBordersx() const
  {
    return m_HistoBorders;
  };
  std::vector<float>& get_HistoBordersx() { return m_HistoBorders; };
  const std::vector<float>& get_HistoBordersy() const
  {
    return m_HistoBordersy;
  };
  std::vector<float>& get_HistoBordersy() { return m_HistoBordersy; };
  const std::vector<float>& get_HistoContents() const
  {
    return m_HistoContents;
  };
  std::vector<float>& get_HistoContents() { return m_HistoContents; };

protected:
  std::vector<float> m_HistoBorders;
  std::vector<float> m_HistoBordersy;
  std::vector<float> m_HistoContents;

private:
  ClassDef(TFCS2DFunctionHistogram, 1)  // TFCS2DFunctionHistogram
};

}  // namespace FastCaloSim::Core

#endif
