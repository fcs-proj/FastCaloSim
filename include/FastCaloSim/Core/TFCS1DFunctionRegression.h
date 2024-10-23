// Copyright (c) 2024 CERN for the benefit of the FastCaloSim project

#ifndef ISF_FASTCALOSIMEVENT_TFCS1DFunctionRegression_h
#define ISF_FASTCALOSIMEVENT_TFCS1DFunctionRegression_h

#include <vector>

#include "FastCaloSim/Core/TFCS1DFunction.h"
#include "TH1.h"
#include "TTree.h"

class TFCS1DFunctionRegression : public TFCS1DFunction
{
public:
  TFCS1DFunctionRegression() {};
  ~TFCS1DFunctionRegression() {};

  using TFCS1DFunction::rnd_to_fct;
  virtual double rnd_to_fct(double rnd) const;
  double regression_value(double uniform) const;
  void set_weights(const std::vector<std::vector<double>>& fWeightMatrix0to1,
                   const std::vector<std::vector<double>>& fWeightMatrix1to2);
  static double sigmoid(double);

private:
  std::vector<std::vector<double>> m_fWeightMatrix0to1;
  std::vector<std::vector<double>> m_fWeightMatrix1to2;

  ClassDef(TFCS1DFunctionRegression, 1)  // TFCS1DFunctionRegression
};

#endif
