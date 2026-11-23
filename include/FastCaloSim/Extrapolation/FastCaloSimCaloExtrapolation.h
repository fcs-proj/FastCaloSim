// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#pragma once

#include <cstddef>
#include <limits>
#include <map>
#include <vector>

#include <FastCaloSim/FastCaloSim_export.h>

#include "FastCaloSim/Core/MLogging.h"

class CaloGeo;
class G4FieldTrack;
class TFCSExtrapolationState;
class TFCSTruthState;

/// @class FastCaloSimCaloExtrapolation
/// @brief Extrapolates a particle, given by the steps of its transport, to the
/// ID-calo boundary and to the entrance, middle and exit of each calorimeter
/// layer.
///
/// Each of these surfaces is approximated by a cylinder around the beam axis.
/// The particle is extrapolated to it along the straight lines between its
/// transport steps.
///
/// The ID-calo boundary and how the barrel layers end depend on the detector:
/// the host sets them with set_settings() before the first extrapolation.
class FASTCALOSIM_EXPORT FastCaloSimCaloExtrapolation
    : public FastCaloSim::MLogging
{
public:
  /// @brief A closed cylinder around the beam axis, centered at the origin.
  /// The dimensions (mm) are used in single precision.
  struct BoundaryCylinder
  {
    double radius;
    double half_length;
  };

  /// @brief The detector-specific part of the extrapolation
  struct Settings
  {
    /// @brief The ID-calo boundary: cylinders, each closer to the beam axis
    /// than the previous one and reaching at least as far in |z|. A particle
    /// crosses the boundary where it reaches one of the cylinders outside the
    /// previous and the next one; of several such points, the closest to the
    /// origin counts. Required.
    std::vector<BoundaryCylinder> id_calo_boundary;

    /// @brief Barrel layers whose cylinders end at the front of another layer:
    /// barrel layer -> other layer. They end at the |z| of the entrance of the
    /// other layer's cells at the largest |eta| on the particle's side (e.g.
    /// the front of the endcap behind the barrel). The cylinders of all other
    /// barrel layers end at the mean |z| of the layer's own ends on both sides.
    std::map<std::size_t, std::size_t> barrel_layer_ends_at_layer;

    /// @brief Particles that cross the ID-calo boundary at this |eta| or
    /// above do not reach the calorimeter and are not extrapolated to the
    /// layers
    double max_abs_eta = std::numeric_limits<double>::infinity();

    /// @brief Failed extrapolations of particles with a pT (MeV) below this
    /// are expected and logged at debug level, the others as warnings
    double min_warning_pt = 0;
  };

  /// @brief Set the detector-specific part of the extrapolation. Throws
  /// std::invalid_argument if the settings are inconsistent, also with the
  /// geometry if it is already set.
  void set_settings(Settings settings);

  auto settings() const -> const Settings& { return m_settings; }

  /// @brief Set the geometry of the calorimeter layers. Throws
  /// std::invalid_argument if it does not fit the settings, if they are
  /// already set.
  void set_geometry(CaloGeo* geo);

  /// @brief Fill result with the positions of the particle on the ID-calo
  /// boundary and in the calorimeter layers. The layers are left empty if the
  /// particle does not reach the calorimeter. Throws std::logic_error if the
  /// settings or the geometry are not set.
  void extrapolate(TFCSExtrapolationState& result,
                   const TFCSTruthState* truth,
                   const std::vector<G4FieldTrack>& caloSteps) const;

private:
  CaloGeo* m_geo = nullptr;
  Settings m_settings;
};
