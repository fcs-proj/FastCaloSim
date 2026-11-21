// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#pragma once

#include <CLHEP/Random/RandGaussZiggurat.h>

namespace FastCaloSim
{
/// Gaussian random number from CLHEP's fast Ziggurat generator. The generator
/// works in single precision: mean and sigma are converted to float and the
/// result back to double, as the parametrizations have always done.
inline auto gauss_ziggurat(CLHEP::HepRandomEngine* engine,
                           double mean,
                           double sigma) -> double
{
  return static_cast<double>(CLHEP::RandGaussZiggurat::shoot(
      engine, static_cast<float>(mean), static_cast<float>(sigma)));
}
}  // namespace FastCaloSim
