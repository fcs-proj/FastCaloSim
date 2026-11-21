// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#include <cmath>

#include "FastCaloSim/Core/TFCS1DFunctionRegression.h"

//=============================================
//======= TFCS1DFunctionRegression =========
//=============================================

/// Network with one input, one hidden layer of sigmoid neurons and one linear
/// output. Row o of m_fWeightMatrix0to1 holds the input weight and the bias of
/// hidden neuron o; m_fWeightMatrix1to2[0] holds the output weights of the
/// hidden neurons followed by the output bias.
auto TFCS1DFunctionRegression::regression_value(double uniform) const -> double
{
  const int n_neurons = static_cast<int>(m_fWeightMatrix0to1.size()) - 1;
  if (n_neurons <= 0)
    return -1;

  const std::vector<double>& output_weights = m_fWeightMatrix1to2[0];
  double result = 0;
  for (int o = 0; o < n_neurons; o++) {
    const std::vector<double>& input_weights = m_fWeightMatrix0to1[o];
    result += output_weights[o]
        * sigmoid(input_weights[0] * uniform + input_weights[1]);
  }
  return result + output_weights[n_neurons];
}

auto TFCS1DFunctionRegression::rnd_to_fct(double rnd) const -> double
{
  return regression_value(rnd);
}

void TFCS1DFunctionRegression::set_weights(
    const std::vector<std::vector<double>>& fWeightMatrix0to1,
    const std::vector<std::vector<double>>& fWeightMatrix1to2)
{
  m_fWeightMatrix0to1 = fWeightMatrix0to1;
  m_fWeightMatrix1to2 = fWeightMatrix1to2;
}

auto TFCS1DFunctionRegression::sigmoid(double x) -> double
{
  return 1 / (1 + std::exp(-x));
}
