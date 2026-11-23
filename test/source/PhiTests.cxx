// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#include <cmath>
#include <numbers>

#include "FastCaloSim/Definitions/Phi.h"

#include <gtest/gtest.h>

namespace
{
constexpr double pi = std::numbers::pi;
}  // namespace

TEST(PhiTests, WrapLeavesAnglesInRangeUnchanged)
{
  for (const double phi : {-pi, -3.0, -1e-300, 0.0, 0.3, 3.0, pi}) {
    EXPECT_EQ(FastCaloSim::wrap_phi(phi), phi);
  }
}

TEST(PhiTests, WrapAnglesOutsideRange)
{
  EXPECT_NEAR(FastCaloSim::wrap_phi(pi + 0.1), -pi + 0.1, 1e-15);
  EXPECT_NEAR(FastCaloSim::wrap_phi(-pi - 0.1), pi - 0.1, 1e-15);
  EXPECT_NEAR(FastCaloSim::wrap_phi(2 * pi + 0.2), 0.2, 1e-15);
  EXPECT_TRUE(std::isnan(FastCaloSim::wrap_phi(std::nan(""))));
}

TEST(PhiTests, DeltaAcrossTheSeam)
{
  // Plain difference away from the seam
  EXPECT_EQ(FastCaloSim::delta_phi(0.5, 0.2), 0.5 - 0.2);
  EXPECT_EQ(FastCaloSim::delta_phi(-3.0, 3.0 - 6.0), -3.0 - (3.0 - 6.0));
  // 3.13 and -3.13 are 2 pi - 6.26 apart across the seam
  EXPECT_NEAR(FastCaloSim::delta_phi(-3.13, 3.13), 2 * pi - 6.26, 1e-14);
  EXPECT_NEAR(FastCaloSim::delta_phi(3.13, -3.13), -(2 * pi - 6.26), 1e-14);
}

TEST(PhiTests, InterpolationAwayFromTheSeamIsThePlainFormula)
{
  for (const double w : {0.0, 0.25, 0.5, 0.7, 1.0}) {
    const double phi1 = -2.9;
    const double phi2 = 0.15;
    EXPECT_EQ(FastCaloSim::interpolate_phi(phi1, phi2, w),
              (1. - w) * phi1 + w * phi2);
  }
}

TEST(PhiTests, InterpolationAcrossTheSeam)
{
  // The midpoint of 3.13 and -3.13 is pi (or -pi), not 0
  EXPECT_NEAR(
      std::abs(FastCaloSim::interpolate_phi(3.13, -3.13, 0.5)), pi, 1e-14);
  // A quarter of the way from 3.13 to -3.13 across the seam
  const double step = 2 * pi - 6.26;
  EXPECT_NEAR(FastCaloSim::interpolate_phi(3.13, -3.13, 0.25),
              3.13 + 0.25 * step,
              1e-14);
  EXPECT_NEAR(FastCaloSim::interpolate_phi(-3.13, 3.13, 0.25),
              -3.13 - 0.25 * step,
              1e-14);
}

TEST(PhiTests, PeriodicDistanceIsTheShorterWayAround)
{
  // Plain distance when the phases are less than half a period apart
  EXPECT_NEAR(FastCaloSim::periodic_distance(0.3, 0.1, 1.0), 0.2, 1e-15);
  EXPECT_NEAR(FastCaloSim::periodic_distance(0.1, 0.3, 1.0), 0.2, 1e-15);
  // Across the end of the period
  EXPECT_NEAR(FastCaloSim::periodic_distance(0.99, 0.01, 1.0), 0.02, 1e-15);
  EXPECT_NEAR(FastCaloSim::periodic_distance(0.01, 0.99, 1.0), 0.02, 1e-15);
  // Half a period apart, the largest distance
  EXPECT_NEAR(FastCaloSim::periodic_distance(0.75, 0.25, 1.0), 0.5, 1e-15);
  // Another period, e.g. the cell size in phi
  const double period = 2 * pi / 256;
  EXPECT_NEAR(
      FastCaloSim::periodic_distance(0.9 * period, 0.1 * period, period),
      0.2 * period,
      1e-15);
  EXPECT_EQ(FastCaloSim::periodic_distance(0.4, 0.4, 1.0), 0.0);
}
