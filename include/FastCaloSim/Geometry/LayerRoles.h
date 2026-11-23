// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#pragma once

#include <algorithm>
#include <cstddef>
#include <functional>
#include <stdexcept>
#include <vector>

namespace FastCaloSim
{

/// @brief The roles that calorimeter layers play in the parametrizations and
/// that the cell geometry does not describe. They depend on the detector, so
/// the host sets them with CaloGeo::set_layer_roles(); a parametrization that
/// needs them throws if they were not set.
struct LayerRoles
{
  /// @brief A layer that the binned showers can take the shower center from
  struct ShowerCenterLayer
  {
    std::size_t layer;
    /// @brief Up to which |eta| of the particle at the layer it is used
    double max_abs_eta;
  };

  /// @brief Strip layers: the electromagnetic layers finely segmented in eta
  /// (ATLAS: EMB1 and EME1). TFCSEnergyAndHitGANV2 places the hits of hadrons
  /// 1 mm apart along the azimuth there instead of 5 mm, and does not limit
  /// them to 10 per voxel and direction.
  std::vector<std::size_t> strip_layers;

  /// @brief The layers the binned showers take the shower center from, in
  /// order: the first one at which the particle's |eta| is not above its
  /// max_abs_eta, the last one (whose max_abs_eta is not used) if there is none
  /// (ATLAS: EMB2 up to |eta| = 1.4, else EME2)
  std::vector<ShowerCenterLayer> shower_center_layers;

  /// @brief Layers in which the upscaling of the binned showers samples the
  /// radius within a sub-bin from a linear distribution instead of a uniform
  /// one (ATLAS: EMB2)
  std::vector<std::size_t> linear_upscaling_layers;

  /// @brief Per layer, the expected energy (MeV) above which
  /// TFCSEnergyRenormalization reports that no cell of the layer received any
  /// of it (ATLAS: 8 times the approximate cell noise, 1500 MeV in EME1 to
  /// EME3). Layers beyond the end are not reported.
  std::vector<double> out_of_calo_report_threshold;

  /// @brief Whether a layer is a strip layer
  auto is_strip_layer(std::size_t layer) const -> bool
  {
    return std::ranges::find(strip_layers, layer) != strip_layers.end();
  }

  /// @brief Whether the binned-shower upscaling samples the radius from a
  /// linear distribution in a layer
  auto is_linear_upscaling_layer(std::size_t layer) const -> bool
  {
    return std::ranges::find(linear_upscaling_layers, layer)
        != linear_upscaling_layers.end();
  }

  /// @brief The layer of the binned-shower center, given the particle's |eta|
  /// at a layer (asked for the candidate layers in order, until one is used)
  auto shower_center_layer(
      const std::function<double(std::size_t)>& abs_eta_at_layer) const
      -> std::size_t
  {
    if (shower_center_layers.empty()) {
      throw std::runtime_error(
          "FastCaloSim::LayerRoles: no shower center layers set");
    }
    for (std::size_t i = 0; i + 1 < shower_center_layers.size(); ++i) {
      const ShowerCenterLayer& candidate = shower_center_layers[i];
      if (abs_eta_at_layer(candidate.layer) > candidate.max_abs_eta) {
        continue;
      }
      return candidate.layer;
    }
    return shower_center_layers.back().layer;
  }
};

}  // namespace FastCaloSim
