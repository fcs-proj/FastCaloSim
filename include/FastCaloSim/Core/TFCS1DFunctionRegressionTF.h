// Copyright (c) 2024 CERN for the benefit of the FastCaloSim project

#ifndef TFCS1DFunctionRegressionTF_h
#define TFCS1DFunctionRegressionTF_h

#include <vector>

#include "FastCaloSim/Core/TFCS1DFunctionRegression.h"
#include "TH1.h"

namespace FastCaloSim::Core
{
class TFCS1DFunctionRegressionTF : public TFCS1DFunctionRegression
{
public:
  TFCS1DFunctionRegressionTF() {};
  TFCS1DFunctionRegressionTF(float, float);
  ~TFCS1DFunctionRegressionTF() {};

  using TFCS1DFunctionRegression::rnd_to_fct;
  virtual double rnd_to_fct(double rnd) const;
  double retransform(double value) const;

private:
  std::vector<std::vector<double>> m_fWeightMatrix0to1;
  std::vector<std::vector<double>> m_fWeightMatrix1to2;
  float m_rangeval;
  float m_startval;

  ClassDef(TFCS1DFunctionRegressionTF, 1)
};

}  // namespace FastCaloSim::Core

#endif
