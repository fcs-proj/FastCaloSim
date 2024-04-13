/*
  Copyright (C) 2002-2019 CERN for the benefit of the ATLAS collaboration
*/

#ifndef ISF_FASTCALOSIMEVENT_TFCS1DFunctionRegressionTF_h
#define ISF_FASTCALOSIMEVENT_TFCS1DFunctionRegressionTF_h

#include <vector>

#include "FastCaloSim/TFCS1DFunctionRegression.h"
#include "TH1.h"

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

#endif
