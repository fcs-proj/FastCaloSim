// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#ifndef ISF_FASTCALOSIMEVENT_TFCS2DFunctionHistogram_h
#define ISF_FASTCALOSIMEVENT_TFCS2DFunctionHistogram_h

#include <vector>

#include "FastCaloSim/Core/TFCS2DFunction.h"

class TH2;

class TFCS2DFunctionHistogram : public TFCS2DFunction
{
public:
  TFCS2DFunctionHistogram(TH2* hist = nullptr)
  {
    if (hist)
      Initialize(hist);
  };
  ~TFCS2DFunctionHistogram() override = default;

  void Initialize(TH2* hist);

  using TFCS2DFunction::rnd_to_fct;
  void rnd_to_fct(float& valuex,
                  float& valuey,
                  float rnd0,
                  float rnd1) const override;

  auto get_HistoBordersx() const -> const std::vector<float>&
  {
    return m_HistoBorders;
  };
  auto get_HistoBordersx() -> std::vector<float>& { return m_HistoBorders; };
  auto get_HistoBordersy() const -> const std::vector<float>&
  {
    return m_HistoBordersy;
  };
  auto get_HistoBordersy() -> std::vector<float>& { return m_HistoBordersy; };
  auto get_HistoContents() const -> const std::vector<float>&
  {
    return m_HistoContents;
  };
  auto get_HistoContents() -> std::vector<float>& { return m_HistoContents; };

protected:
  std::vector<float> m_HistoBorders;
  std::vector<float> m_HistoBordersy;
  std::vector<float> m_HistoContents;

private:
  ClassDef(TFCS2DFunctionHistogram, 1)  // TFCS2DFunctionHistogram
};

#endif
