// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#pragma once

#include <cmath>
#include <numbers>

namespace FastCaloSim
{
/// The azimuth phi in [-pi, pi]; unchanged if it is already in that range
inline auto wrap_phi(double phi) -> double
{
  if (std::abs(phi) <= std::numbers::pi) {
    return phi;
  }
  return std::remainder(phi, 2 * std::numbers::pi);
}

/// The azimuthal difference phi1 - phi2 in [-pi, pi]: the plain difference
/// unless the two angles lie more than pi apart, i.e. on different sides of
/// the +-pi seam
inline auto delta_phi(double phi1, double phi2) -> double
{
  return wrap_phi(phi1 - phi2);
}

/// The azimuth (1 - w) * phi1 + w * phi2 along the shorter arc between phi1
/// and phi2, in [-pi, pi]: the plain interpolation unless the two angles lie
/// on different sides of the +-pi seam
inline auto interpolate_phi(double phi1, double phi2, double w) -> double
{
  if (std::abs(phi2 - phi1) <= std::numbers::pi) {
    return (1. - w) * phi1 + w * phi2;
  }
  return wrap_phi((1. - w) * phi1 + w * (phi1 + wrap_phi(phi2 - phi1)));
}

/// The distance between two phases a and b of a quantity with the given
/// period, in [0, period / 2]: the shorter way around, e.g. 0.02 for phases
/// 0.99 and 0.01 of a period of 1
inline auto periodic_distance(double a, double b, double period) -> double
{
  return std::abs(std::remainder(a - b, period));
}
}  // namespace FastCaloSim
