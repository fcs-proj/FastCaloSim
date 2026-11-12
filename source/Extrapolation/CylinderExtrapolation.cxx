// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#include <cmath>

#include "CylinderExtrapolation.h"

namespace FastCaloSim::Extrapolation
{

namespace
{

/// Distance between the point and the segment from a to b
auto distanceToSegment(const Vector3D& point,
                       const Vector3D& a,
                       const Vector3D& b) -> double
{
  const Vector3D ab = b - a;
  const double c1 = (point - a).dot(ab);
  if (c1 <= 0) {
    return (a - point).mag();
  }
  const double c2 = ab.dot(ab);
  if (c2 <= c1) {
    return (b - point).mag();
  }
  return (a + (c1 / c2) * ab - point).mag();
}

/// Whether the point is on the segment from a to b, within 1 um
auto isOnSegment(const Vector3D& point, const Vector3D& a, const Vector3D& b)
    -> bool
{
  return distanceToSegment(point, a, b) < 0.001;
}

/// The point at parameter t on the line from a (t = 0) to b (t = 1), clamped
/// to the segment between them
auto pointOnSegment(const Vector3D& a, const Vector3D& b, double t) -> Vector3D
{
  if (t <= 0) {
    return a;
  }
  if (t >= 1) {
    return b;
  }
  return a + t * (b - a);
}

/// The point of the line through posA along dirA closest to the line through
/// posB along dirB; posA if the lines are parallel
auto closestPointOfLines(const Vector3D& posA,
                         const Vector3D& dirA,
                         const Vector3D& posB,
                         const Vector3D& dirB) -> Vector3D
{
  const Vector3D ea = dirA.unit();
  const Vector3D eb = dirB.unit();
  const Vector3D ab = posB - posA;
  const double eaDotEb = ea.dot(eb);
  const double denominator = 1 - eaDotEb * eaDotEb;
  if (std::abs(denominator) > 10e-7) {
    return posA + ((ab.dot(ea) - ab.dot(eb) * eaDotEb) / denominator) * ea;
  }
  return posA;
}

/// The intersections of the infinite line through a and b, projected onto the
/// xy plane, with the circle of the given radius around the z axis. For lines
/// parallel to the endcaps: the points take the z of a (first) and b (second).
auto circleIntersections(float radius, const Vector3D& a, const Vector3D& b)
    -> Intersections
{
  const double dx = b.x() - a.x();
  const double dy = b.y() - a.y();
  const double qa = dx * dx + dy * dy;
  const double qb = 2 * (dx * a.x() + dy * a.y());
  const double qc = a.x() * a.x() + a.y() * a.y() - radius * radius;
  const double discriminant = qb * qb - 4 * qa * qc;

  const auto pointAt = [&](double t, double z)
  { return Vector3D(a.x() + t * dx, a.y() + t * dy, z); };

  Intersections result;
  if (qa <= 0.0000001 || discriminant < 0) {
    return result;
  }
  if (std::abs(discriminant) < 0.00001) {
    // The line touches the circle
    result.points[0] = pointAt(-qb / (2 * qa), a.z());
    result.count = 1;
    return result;
  }
  result.points = {pointAt((-qb + std::sqrt(discriminant)) / (2 * qa), a.z()),
                   pointAt((-qb - std::sqrt(discriminant)) / (2 * qa), b.z())};
  result.count = 2;
  return result;
}

/// The point between the bounds a and b closest to the surface, found by
/// walking from both bounds towards each other in steps of 10 um until the
/// distance to the surface grows
auto iterativeClosestPoint(const Cylinder& cylinder,
                           const Vector3D& boundA,
                           const Vector3D& boundB) -> Vector3D
{
  constexpr double stepSize = 0.01;  // mm

  const Vector3D boundDir = boundB - boundA;
  const double distBounds = boundDir.mag();
  // Nothing to search between bounds this close
  if (distBounds <= 4 * stepSize) {
    return boundA + 0.5 * (boundB - boundA);
  }
  const Vector3D unitDir = boundDir / distBounds;

  double minDistA = cylinder.distanceTo(boundA);
  double minDistB = cylinder.distanceTo(boundB);
  // Without a minimum in between, the bound closer to the surface
  Vector3D resultA = minDistA < minDistB ? boundA : boundB;
  Vector3D resultB = resultA;

  const auto nHalfDivisions =
      static_cast<unsigned int>((distBounds / stepSize) / 2);
  for (unsigned int step = 0; step < nHalfDivisions; step++) {
    const double previous = (step * stepSize) / distBounds;

    const Vector3D tmpA = boundA + (step + 1) * stepSize * unitDir;
    const double tmpDistA = cylinder.distanceTo(tmpA);
    if (minDistA < tmpDistA) {
      resultA = boundA + previous * boundDir;
      resultB = tmpA;
      break;
    }
    minDistA = tmpDistA;

    const Vector3D tmpB = boundB - (step + 1) * stepSize * unitDir;
    const double tmpDistB = cylinder.distanceTo(tmpB);
    if (minDistB < tmpDistB) {
      resultB = boundB - previous * boundDir;
      resultA = tmpB;
      break;
    }
    minDistB = tmpDistB;
  }

  // The middle of the step in which the distance started to grow
  return resultA + 0.5 * (resultB - resultA);
}

/// The point of the segment from a to b closest to the surface (approximately)
auto closestPointToSurface(const Cylinder& cylinder,
                           const Vector3D& a,
                           const Vector3D& b) -> Vector3D
{
  const float radius = cylinder.radius;
  const float halfLength = cylinder.halfLength;
  const Vector3D ab = b - a;
  const Vector3D projA = cylinder.project(a);
  const Vector3D projB = cylinder.project(b);

  // Both ends project onto the same point of the surface
  if ((projB - projA).mag() < 0.0001) {
    return a;
  }

  const bool beyondEndcaps = a.z() > halfLength || a.z() < -halfLength
      || b.z() > halfLength || b.z() < -halfLength;
  if (!beyondEndcaps) {
    // The point closest to the z axis if on the segment, else the end closer
    // to the surface
    const Vector3D closestToAxis =
        closestPointOfLines(a, ab.unit(), Vector3D(0, 0, 0), Vector3D(0, 0, 1));
    if (isOnSegment(closestToAxis, a, b)) {
      return closestToAxis;
    }
    return (projA - a).mag() < (projB - b).mag() ? a : b;
  }

  // The endcap on the side of the segment
  const double endcapZ = a.z() + b.z() >= 0 ? halfLength : -halfLength;
  const Vector3D endcapCenter(0, 0, endcapZ);
  const bool parallelToEndcaps = std::abs(a.z() - b.z()) < 0.00001;
  if (!parallelToEndcaps) {
    // Search between the point closest to the endcap center and the crossing
    // of the endcap plane
    const Vector3D closestToCenter =
        pointOnSegment(a, b, (endcapCenter - a).dot(ab) / ab.dot(ab));
    const Vector3D planeCrossing =
        pointOnSegment(a, b, (endcapZ - a.z()) / ab.z());
    return iterativeClosestPoint(cylinder, closestToCenter, planeCrossing);
  }

  // The segment is parallel to the endcaps and beyond them, so all its points
  // within the radius are equally close to the endcap disc: take the middle of
  // the part within the radius
  const Intersections circle = circleIntersections(radius, a, b);
  if (circle.count == 2) {
    const auto& [p, q] = circle.points;
    const bool pOnSegment = isOnSegment(p, a, b);
    const bool qOnSegment = isOnSegment(q, a, b);
    if (pOnSegment && qOnSegment) {
      return p + 0.5 * (q - p);
    }
    const Vector3D& inner = a.perp() <= radius ? a : b;
    if (pOnSegment) {
      return p + 0.5 * (inner - p);
    }
    if (qOnSegment) {
      return q + 0.5 * (inner - q);
    }
    return a + 0.5 * ab;
  }
  if (circle.count == 1) {
    // The line touches the circle: there, if on the segment, else at the end
    // closer to it
    const Vector3D& touch = circle.points[0];
    if (isOnSegment(touch, a, b)) {
      return touch;
    }
    return (a - touch).mag() < (b - touch).mag() ? a : b;
  }
  // The line misses the circle: the point closest to the endcap center if on
  // the segment, else the end closer to it
  const Vector3D closestToCenter =
      a + ((endcapCenter - a).dot(ab) / ab.dot(ab)) * ab;
  if (isOnSegment(closestToCenter, a, b)) {
    return closestToCenter;
  }
  return (a - closestToCenter).mag() < (b - closestToCenter).mag() ? a : b;
}

/// The momentum direction of the track point closer to the surface
auto closerDirection(const Cylinder& cylinder,
                     const TrackPoint& start,
                     const TrackPoint& end) -> const Vector3D&
{
  return cylinder.distanceTo(start.position) < cylinder.distanceTo(end.position)
      ? start.direction
      : end.direction;
}

/// Of the two intersections of the line through the segment, the one the
/// particle crosses
auto selectIntersection(const Vector3D& start,
                        const Vector3D& end,
                        Location startLocation,
                        Location endLocation,
                        const Intersections& intersections) -> const Vector3D&
{
  const auto& [a, b] = intersections.points;
  const bool startInside = startLocation == Location::Inside;
  const bool endInside = endLocation == Location::Inside;
  if (startInside != endInside) {
    // The segment crosses the surface: the intersection closer to it
    return distanceToSegment(a, start, end) > distanceToSegment(b, start, end)
        ? b
        : a;
  }
  if (startInside) {
    // Both inside: the intersection further ahead
    const Vector3D direction = end - start;
    return (a - end).dot(direction) < (b - end).dot(direction) ? b : a;
  }
  // Both outside: the intersection closer to the end
  return (end - a).mag() > (end - b).mag() ? b : a;
}

/// Where the track crosses the surface: the first segment that crosses it or
/// whose line intersects it on the segment, with a few exceptions for
/// intersections behind the segment, else the intersection of the line of the
/// last segment
auto extrapolateWithIntersection(const std::vector<TrackPoint>& track,
                                 const Cylinder& cylinder)
    -> std::optional<TrackPoint>
{
  // Number of segments whose intersections with the surface were computed
  unsigned int nExtrapolations = 0;
  for (std::size_t i = 1; i < track.size(); ++i) {
    const TrackPoint& start = track[i - 1];
    const TrackPoint& end = track[i];
    const Vector3D& p1 = start.position;
    const Vector3D& p2 = end.position;

    const Location location1 = cylinder.locate(p1);
    const Location location2 = cylinder.locate(p2);
    if (location1 == Location::OnSurface) {
      return start;
    }
    if (location2 == Location::OnSurface) {
      return end;
    }

    // Skip segments too short for a reliable line
    if ((p2 - p1).mag() < 0.01) {
      continue;
    }

    const Intersections intersections = cylinder.intersections(p1, p2);
    nExtrapolations++;
    if (intersections.count == 0) {
      continue;
    }
    const Vector3D& intersection = intersections.count == 1
        ? intersections.points[0]
        : selectIntersection(p1, p2, location1, location2, intersections);

    const bool isForward = (intersection - p1).dot(p2 - p1) >= 0;
    const bool crossesSurface =
        (location1 == Location::Inside) != (location2 == Location::Inside);

    // Only the first segment may intersect behind it, unless it crosses
    if (nExtrapolations > 1 && !isForward && !crossesSurface) {
      continue;
    }

    const bool bothOutside =
        location1 == Location::Outside && location2 == Location::Outside;
    const bool isLast = i == track.size() - 1;
    if (crossesSurface || isOnSegment(intersection, p1, p2)
        || (bothOutside && !isForward && nExtrapolations == 1) || isLast)
    {
      return TrackPoint {intersection, closerDirection(cylinder, start, end)};
    }
  }
  return std::nullopt;
}

/// The point of the surface closest to the track, if within 100 m
auto extrapolateWithClosestApproach(const std::vector<TrackPoint>& track,
                                    const Cylinder& cylinder)
    -> std::optional<TrackPoint>
{
  std::optional<TrackPoint> result;
  double minDistance = 100000;  // mm
  for (std::size_t i = 1; i < track.size(); ++i) {
    const TrackPoint& start = track[i - 1];
    const TrackPoint& end = track[i];
    const Vector3D closest =
        closestPointToSurface(cylinder, start.position, end.position);
    const Vector3D onSurface = cylinder.project(closest);
    const double distance = (closest - onSurface).mag();
    if (distance < minDistance) {
      minDistance = distance;
      result = TrackPoint {onSurface, closerDirection(cylinder, start, end)};
    }
  }
  return result;
}

}  // namespace

auto extrapolateToCylinder(const std::vector<TrackPoint>& track,
                           const Cylinder& cylinder)
    -> std::optional<TrackPoint>
{
  if (track.size() == 1) {
    return TrackPoint {cylinder.project(track[0].position), track[0].direction};
  }
  if (auto crossing = extrapolateWithIntersection(track, cylinder)) {
    return crossing;
  }
  return extrapolateWithClosestApproach(track, cylinder);
}

}  // namespace FastCaloSim::Extrapolation
