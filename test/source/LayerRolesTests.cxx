// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

#include "FastCaloSim/Geometry/LayerRoles.h"

#include <gtest/gtest.h>

#include "FastCaloSim/Geometry/CaloGeo.h"

namespace
{
// Two candidate layers of the shower center, as in ATLAS
auto two_shower_center_layers() -> FastCaloSim::LayerRoles
{
  FastCaloSim::LayerRoles roles;
  roles.shower_center_layers = {
      {.layer = 2, .max_abs_eta = 1.4},
      {.layer = 6, .max_abs_eta = 0},  // not used: the last one
  };
  return roles;
}
}  // namespace

// Without roles set by the host, asking for them throws
TEST(LayerRolesTests, CaloGeoWithoutRoles)
{
  const CaloGeo geo;
  EXPECT_FALSE(geo.has_layer_roles());
  EXPECT_THROW((void)geo.layer_roles(), std::runtime_error);
}

TEST(LayerRolesTests, CaloGeoWithRoles)
{
  CaloGeo geo;
  FastCaloSim::LayerRoles roles;
  roles.strip_layers = {1, 5};
  geo.set_layer_roles(roles);
  ASSERT_TRUE(geo.has_layer_roles());
  EXPECT_EQ(geo.layer_roles().strip_layers, (std::vector<std::size_t> {1, 5}));
}

TEST(LayerRolesTests, StripAndLinearUpscalingLayers)
{
  FastCaloSim::LayerRoles roles;
  roles.strip_layers = {1, 5};
  roles.linear_upscaling_layers = {2};
  EXPECT_TRUE(roles.is_strip_layer(1));
  EXPECT_TRUE(roles.is_strip_layer(5));
  EXPECT_FALSE(roles.is_strip_layer(2));
  EXPECT_TRUE(roles.is_linear_upscaling_layer(2));
  EXPECT_FALSE(roles.is_linear_upscaling_layer(6));
}

// The first candidate layer at which |eta| is not above its limit, else the
// last one
TEST(LayerRolesTests, ShowerCenterLayer)
{
  const FastCaloSim::LayerRoles roles = two_shower_center_layers();
  const auto at = [&roles](double abs_eta_at_layer_2) -> std::size_t
  {
    return roles.shower_center_layer(
        [abs_eta_at_layer_2](std::size_t layer) -> double
        {
          EXPECT_EQ(layer, 2U);  // the last candidate is never asked for
          return abs_eta_at_layer_2;
        });
  };
  EXPECT_EQ(at(0.0), 2U);
  EXPECT_EQ(at(1.4), 2U);
  EXPECT_EQ(at(1.41), 6U);
  EXPECT_EQ(at(999.0), 6U);
  // NaN is not above the limit
  EXPECT_EQ(at(std::numeric_limits<double>::quiet_NaN()), 2U);
}

TEST(LayerRolesTests, NoShowerCenterLayers)
{
  const FastCaloSim::LayerRoles roles;
  EXPECT_THROW(
      (void)roles.shower_center_layer([](std::size_t) -> double { return 0; }),
      std::runtime_error);
}
