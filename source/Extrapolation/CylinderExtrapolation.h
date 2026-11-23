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

/// @brief Where a track reaches the surface of a cylinder.
///
/// The track is given by its points after each transport step and is a
/// straight line between consecutive points. The result is a track point on
/// the surface, or where the track, or its straight continuation before its
/// first or after its last point, crosses the surface. Its direction is the one
/// of the track point closer to the surface. There is no result for a track
/// with fewer than two points or whose continuation misses the surface.
auto crossCylinder(const std::vector<TrackPoint>& track,
                   const Cylinder& cylinder) -> std::optional<TrackPoint>;

/// @brief Extrapolates a track to the surface of a cylinder.
///
/// The result is where the track reaches the surface (see crossCylinder) or,
/// if it does not, the point of the surface closest to the track. Its direction
/// is the one of the track point closer to the surface. There is no result
/// only for a track without points.
auto extrapolateToCylinder(const std::vector<TrackPoint>& track,
                           const Cylinder& cylinder)
    -> std::optional<TrackPoint>;

}  // namespace FastCaloSim::Extrapolation
