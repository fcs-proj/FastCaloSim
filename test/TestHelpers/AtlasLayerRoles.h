// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#pragma once

#include <array>
#include <cstddef>
#include <initializer_list>
#include <limits>
#include <numbers>
#include <vector>

#include "FastCaloSim/Geometry/LayerRoles.h"

namespace TestHelpers
{

/// The roles of the ATLAS calorimeter layers in the parametrizations, as the
/// ATLAS simulation sets them
inline auto atlas_layer_roles() -> FastCaloSim::LayerRoles
{
  // The ATLAS calorimeter layers (CaloCell_ID::CaloSample)
  enum : std::size_t
  {
    PreSamplerB,
    EMB1,
    EMB2,
    EMB3,
    PreSamplerE,
    EME1,
    EME2,
    EME3,
    HEC0,
    HEC1,
    HEC2,
    HEC3,
    TileBar0,
    TileBar1,
    TileBar2,
    TileGap1,
    TileGap2,
    TileGap3,
    TileExt0,
    TileExt1,
    TileExt2,
    FCAL0,
    FCAL1,
    FCAL2,
    nLayers
  };

  // Approximate cell noise per layer (MeV)
  constexpr std::array<double, nLayers> approxLayerNoise {
      150.0, 40.0,  80.0,  40.0,  150.0, 40.0,  80.0,  50.0,
      400.0, 400.0, 400.0, 400.0, 300.0, 150.0, 40.0,  150.0,
      40.0,  400.0, 300.0, 150.0, 40.0,  400.0, 400.0, 400.0};

  FastCaloSim::LayerRoles roles;
  roles.strip_layers = {EMB1, EME1};
  roles.shower_center_layers = {
      {.layer = EMB2, .max_abs_eta = 1.4},
      {.layer = EME2, .max_abs_eta = std::numeric_limits<double>::infinity()}};
  roles.linear_upscaling_layers = {EMB2};
  roles.out_of_calo_report_threshold.resize(nLayers);
  for (std::size_t layer = 0; layer < nLayers; ++layer) {
    roles.out_of_calo_report_threshold[layer] = 8 * approxLayerNoise[layer];
  }
  for (const std::size_t layer : {EME1, EME2, EME3}) {
    roles.out_of_calo_report_threshold[layer] = 1500;
  }
  // The accordion absorbers: 1024 in phi in the barrel, 768 in the outer and
  // 256 in the inner wheel of the endcap
  using Period = FastCaloSim::LayerRoles::PhiModulationPeriod;
  constexpr double inf = std::numeric_limits<double>::infinity();
  const std::vector<Period> barrel {
      {.max_abs_eta = inf, .period = 2 * std::numbers::pi / 1024}};
  const std::vector<Period> endcap {
      {.max_abs_eta = 2.5, .period = 2 * std::numbers::pi / 768},
      {.max_abs_eta = inf, .period = 2 * std::numbers::pi / 256}};
  roles.phi_modulation_periods.resize(nLayers);
  for (std::size_t layer = 0; layer < nLayers; ++layer) {
    roles.phi_modulation_periods[layer] = layer <= EMB3 ? barrel : endcap;
  }
  return roles;
}

}  // namespace TestHelpers
