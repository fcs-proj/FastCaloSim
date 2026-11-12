// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#pragma once

#include <optional>
#include <vector>

#include "Cylinder.h"

namespace FastCaloSim::Extrapolation
{

/// @brief Position and momentum direction of a particle
struct TrackPoint
{
  Vector3D position;
  Vector3D direction;
};

/// @brief Extrapolates a track to the surface of a cylinder.
///
/// The track is given by its points after each transport step and is a
/// straight line between consecutive points. The result is where the track
/// crosses the surface or, if it does not, the point of the surface closest to
/// the track. Its direction is the one of the track point closer to the
/// surface. There is no result only for a track without points.
auto extrapolateToCylinder(const std::vector<TrackPoint>& track,
                           const Cylinder& cylinder)
    -> std::optional<TrackPoint>;

}  // namespace FastCaloSim::Extrapolation
