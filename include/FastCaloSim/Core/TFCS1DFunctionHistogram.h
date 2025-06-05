// Copyright (c) 2024 CERN for the benefit of the FastCaloSim project

#ifndef TFCS1DFunctionHistogram_h
#define TFCS1DFunctionHistogram_h

#include <vector>

#include "FastCaloSim/Core/TFCS1DFunction.h"
#include "TH1.h"

namespace FastCaloSim::Core
{
class TFCS1DFunctionHistogram : public TFCS1DFunction
{
public:
  TFCS1DFunctionHistogram() {};
  TFCS1DFunctionHistogram(TH1* hist, double);

  void Initialize(TH1* hist, double);

  using TFCS1DFunction::rnd_to_fct;
  virtual double rnd_to_fct(double rnd) const;
  TH1* vector_to_histo();
  double get_inverse(double rnd) const;
  static double linear(double y1, double y2, double x1, double x2, double x);
  static double non_linear(
      double y1, double y2, double x1, double x2, double x);

  static double get_maxdev(TH1*, TH1D*);
  void smart_rebin_loop(TH1* hist, double);
  static double get_change(TH1*);
  static TH1D* smart_rebin(TH1D*);
  static double* histo_to_array(TH1*);
  static double sample_from_histo(TH1* hist, double);
  double sample_from_histovalues(double);

  std::vector<float> get_HistoBorders() { return m_HistoBorders; };
  std::vector<float> get_HistoContents() { return m_HistoContents; };

protected:
  std::vector<float> m_HistoBorders;
  std::vector<float> m_HistoContents;

  ClassDef(TFCS1DFunctionHistogram, 1)  // TFCS1DFunctionHistogram
};
}  // namespace FastCaloSim::Core

#endif
