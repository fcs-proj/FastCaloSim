// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#include <algorithm>
#include <cmath>
#include <optional>

#include "Cylinder.h"

namespace FastCaloSim::Extrapolation
{

namespace
{

/// The intersections of the infinite line through a and b with the cover of
/// the cylinder (without the endcaps)
auto coverIntersections(const Cylinder& cylinder,
                        const Vector3D& a,
                        const Vector3D& b) -> Intersections
{
  // Find where the line, projected onto the xy plane, crosses the circle of
  // the cylinder radius, from the distance of the line to the z axis
  const Vector3D aXY(a.x(), a.y(), 0);
  const Vector3D bXY(b.x(), b.y(), 0);
  const Vector3D diffXY = aXY - bXY;
  const double norm2 = diffXY.dot(diffXY);
  // A line parallel to the z axis does not cross the cover in a point
  if (norm2 == 0) {
    return {};
  }
  const double t = aXY.dot(diffXY) / norm2;
  // Squared distance to the z axis, which rounding can make negative
  const double d2 = std::max(aXY.dot(aXY) - t * t * norm2, 0.0);

  const float radius = cylinder.radius;
  if (d2 > radius * radius) {
    return {};
  }
  const double k = std::sqrt((radius * radius - d2) / norm2);

  // Keep the crossings within the z range of the cover
  Intersections result;
  for (const double s : {t + k, t - k}) {
    const Vector3D point = a + s * (b - a);
    if (point.z() <= cylinder.halfLength && point.z() >= -cylinder.halfLength) {
      result.points[result.count++] = point;
    }
  }
  return result;
}

/// The intersection of the infinite line through a and b with the endcap at
/// z = endcapZ, if there is one
auto endcapIntersection(const Cylinder& cylinder,
                        const Vector3D& a,
                        const Vector3D& b,
                        double endcapZ) -> std::optional<Vector3D>
{
  const Vector3D ab = b - a;
  // Lines (almost) parallel to the endcap do not cross it
  if (std::abs(ab.z()) <= 1e-6) {
    return std::nullopt;
  }
  const Vector3D center(0, 0, endcapZ);
  const Vector3D intersection = b + ((center - b).z() / ab.z()) * ab;
  if ((intersection - center).mag() <= cylinder.radius) {
    return intersection;
  }
  return std::nullopt;
}

}  // namespace

auto Cylinder::locate(const Vector3D& point) const -> Location
{
  // Larger tolerances lead to extrapolation failures around |eta| = 4
  constexpr float tolerance = 1;  // mm

  const double perp = point.perp();
  const double z = point.z();
  const bool withinZ = z < halfLength && z > -halfLength;
  const double distanceToEndcap =
      z > 0 ? std::abs(z - halfLength) : std::abs(z + halfLength);

  const bool onEndcap =
      perp <= radius + tolerance && distanceToEndcap < tolerance;
  const bool onCover = std::abs(perp - radius) < tolerance && withinZ;
  if (onEndcap || onCover) {
    return Location::OnSurface;
  }
  if (withinZ && perp < radius) {
    return Location::Inside;
  }
  return Location::Outside;
}

auto Cylinder::project(const Vector3D& point) const -> Vector3D
{
  const double z = point.z();
  if (z >= halfLength || z <= -halfLength) {
    // Beyond an endcap: onto the endcap disc, at the edge if outside of it
    const Vector3D center(0, 0, z >= halfLength ? halfLength : -halfLength);
    const Vector3D inEndcapPlane(point.x(), point.y(), center.z());
    return center
        + std::min(point.perp(), static_cast<double>(radius))
        * (inEndcapPlane - center).unit();
  }
  // Within the z range: radially onto the cover
  const Vector3D onAxis(0, 0, z);
  return onAxis + radius * (point - onAxis).unit();
}

auto Cylinder::intersections(const Vector3D& a, const Vector3D& b) const
    -> Intersections
{
  Intersections result = coverIntersections(*this, a, b);
  if (result.count == 2) {
    return result;
  }

  const auto positive = endcapIntersection(*this, a, b, halfLength);
  const auto negative = endcapIntersection(*this, a, b, -halfLength);
  if (positive && negative) {
    if (result.count == 1) {
      // The cover intersection is on the edge to one of the endcaps, which
      // the line then crosses at the same point: take the farther one
      result.points[1] = (*positive - result.points[0]).mag()
              > (*negative - result.points[0]).mag()
          ? *positive
          : *negative;
    } else {
      result.points = {*positive, *negative};
    }
    result.count = 2;
  } else if (positive || negative) {
    result.points[result.count++] = positive ? *positive : *negative;
  }
  return result;
}

}  // namespace FastCaloSim::Extrapolation
