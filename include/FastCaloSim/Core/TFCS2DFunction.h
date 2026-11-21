// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#pragma once

#include <vector>

#include "FastCaloSim/Core/TFCSFunction.h"

class TH2;

class TFCS2DFunction : public TFCSFunction
{
public:
  TFCS2DFunction() = default;
  ~TFCS2DFunction() override = default;
  TFCS2DFunction(const TFCS2DFunction&) = default;
  TFCS2DFunction(TFCS2DFunction&&) = default;
  auto operator=(const TFCS2DFunction&) -> TFCS2DFunction& = default;
  auto operator=(TFCS2DFunction&&) -> TFCS2DFunction& = default;

  auto ndim() const -> int override { return 2; }

  virtual void rnd_to_fct(double& valuex,
                          double& valuey,
                          double rnd0,
                          double rnd1) const = 0;
  void rnd_to_fct(float* value, const float* rnd) const override;

  static auto CheckAndIntegrate2DHistogram(const TH2* hist,
                                           std::vector<double>& integral_vec,
                                           int& first,
                                           int& last) -> double;

private:
  ClassDefOverride(TFCS2DFunction, 1)  // TFCS2DFunction
};
