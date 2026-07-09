// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#pragma once

#include <vector>

#include "FastCaloSim/Core/TFCS1DFunctionRegression.h"
#include "TH1.h"

class TFCS1DFunctionRegressionTF : public TFCS1DFunctionRegression
{
public:
  TFCS1DFunctionRegressionTF() = default;
  TFCS1DFunctionRegressionTF(float, float);
  ~TFCS1DFunctionRegressionTF() override = default;

  using TFCS1DFunctionRegression::rnd_to_fct;
  auto rnd_to_fct(double rnd) const -> double override;
  auto retransform(double value) const -> double;

private:
  std::vector<std::vector<double>> m_fWeightMatrix0to1;
  std::vector<std::vector<double>> m_fWeightMatrix1to2;
  float m_rangeval {};
  float m_startval {};

  ClassDefOverride(TFCS1DFunctionRegressionTF, 1)
};
