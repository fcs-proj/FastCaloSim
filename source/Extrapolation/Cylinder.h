// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#pragma once

#include <array>

#include <CLHEP/Vector/ThreeVector.h>

namespace FastCaloSim::Extrapolation
{

using Vector3D = CLHEP::Hep3Vector;

/// @brief Where a point lies with respect to a cylinder
enum class Location
{
  Inside,
  Outside,
  OnSurface
};

/// @brief The points where a line crosses a surface
struct Intersections
{
  std::array<Vector3D, 2> points;
  int count = 0;
};

/// @brief A closed cylinder around the z axis, centered at the origin
struct Cylinder
{
  /// The dimensions are kept in single precision, as the extrapolation has
  /// always done: double precision would move the extrapolated positions.
  constexpr Cylinder(double radius, double halfLength)
      : radius(static_cast<float>(radius))
      , halfLength(static_cast<float>(halfLength))
  {
  }

  float radius;
  float halfLength;

  /// @brief Where the point lies; points within 1 mm of the surface are on it
  auto locate(const Vector3D& point) const -> Location;

  /// @brief Projects the point onto the surface: onto the endcap disc for
  /// points beyond an endcap, radially onto the cover for all others
  auto project(const Vector3D& point) const -> Vector3D;

  /// @brief Distance between the point and its projection onto the surface
  auto distanceTo(const Vector3D& point) const -> double
  {
    return (point - project(point)).mag();
  }

  /// @brief The intersections of the infinite line through a and b with the
  /// surface (cover and endcaps)
  auto intersections(const Vector3D& a, const Vector3D& b) const
      -> Intersections;
};

}  // namespace FastCaloSim::Extrapolation
