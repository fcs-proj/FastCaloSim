// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#include <array>
#include <cmath>
#include <optional>
#include <string>

#include "FastCaloSim/Extrapolation/FastCaloSimCaloExtrapolation.h"

#include <fmt/format.h>

#include "CylinderExtrapolation.h"
#include "FastCaloSim/Core/TFCSExtrapolationState.h"
#include "FastCaloSim/Core/TFCSTruthState.h"
#include "FastCaloSim/Geometry/CaloGeo.h"
#include "G4FieldTrack.hh"

using FastCaloSim::Extrapolation::Cylinder;
using FastCaloSim::Extrapolation::extrapolateToCylinder;
using FastCaloSim::Extrapolation::TrackPoint;
using FastCaloSim::Extrapolation::Vector3D;

namespace
{

/// The ID-calo boundary: three cylinders, each closer to the beam axis and
/// reaching at least as far in |z| as the previous one
constexpr std::array<Cylinder, 3> kIDCaloBoundary {
    {{1148.0, 3550.0}, {120.0, 4587.0}, {41.0, 4587.0}}};

/// Failed extrapolations of particles above this pT (MeV) are warnings, those
/// of others are expected and logged at debug level
constexpr double kWarningPtLimit = 500;

/// Where the track leaves the ID: of the points where it reaches the boundary
/// cylinders, the closest to the origin that is on the boundary, i.e. neither
/// within the z range of the previous cylinder nor the radius of the next one
auto extrapolateToIDCaloBoundary(const std::vector<TrackPoint>& track)
    -> std::optional<TrackPoint>
{
  constexpr double tolerance = 0.001;  // mm

  std::optional<TrackPoint> boundary;
  for (std::size_t i = 0; i < kIDCaloBoundary.size(); ++i) {
    const auto point = extrapolateToCylinder(track, kIDCaloBoundary[i]);
    if (!point) {
      continue;
    }
    const Vector3D& position = point->position;
    const bool withinPrevious = i > 0
        && std::abs(position.z())
            < kIDCaloBoundary[i - 1].halfLength + tolerance;
    const bool withinNext = i + 1 < kIDCaloBoundary.size()
        && position.perp() < kIDCaloBoundary[i + 1].radius + tolerance;
    if (withinPrevious || withinNext) {
      continue;
    }
    if (!boundary || position.mag() <= boundary->position.mag()) {
      boundary = point;
    }
  }
  return boundary;
}

/// The three cylinders of a layer, indexed by Cell::SubPos
template<typename CylinderOfSubPos>
auto perSubPos(CylinderOfSubPos&& cylinder) -> std::array<Cylinder, 3>
{
  static_assert(Cell::SubPos::MID == 0 && Cell::SubPos::ENT == 1
                && Cell::SubPos::EXT == 2);
  return {cylinder(Cell::SubPos::MID),
          cylinder(Cell::SubPos::ENT),
          cylinder(Cell::SubPos::EXT)};
}

/// The cylinders through the middle, entrance and exit of the calorimeter
/// layers, for a particle that crossed the ID-calo boundary at a position.
/// They go through the cell at that position: barrel layers are cylinder
/// covers, closed at the end of the layer; the other layers are endcaps, out to
/// the edge of the layer closest to the barrel.
class LayerCylinders
{
public:
  LayerCylinders(const CaloGeo& geo, const Vector3D& boundary)
      : m_geo(geo)
      , m_boundary {boundary.x(),
                    boundary.y(),
                    boundary.z(),
                    boundary.eta(),
                    boundary.phi(),
                    boundary.perp()}
      , m_farPositive {0, 0, 0, 1000, boundary.phi(), 0}
      , m_farNegative {0, 0, 0, -1000, boundary.phi(), 0}
      , m_positiveSide(boundary.eta() > 0)
  {
  }

  /// The cylinders of the layer, indexed by Cell::SubPos
  auto of(unsigned int layer) -> std::array<Cylinder, 3>
  {
    const Cell cell = m_geo.get_cell(layer, m_boundary);

    if (m_geo.is_barrel(layer)) {
      if (layer < 4) {
        // EMB: closed at the front of EME1 on the particle's side
        const double halfLength = embHalfLength();
        return perSubPos(
            [&](Cell::SubPos subpos)
            { return Cylinder(std::abs(cell.r(subpos)), halfLength); });
      }
      // Tile: closed at the mean |z| of the layer's ends on both sides
      const Cell positiveEnd = m_geo.get_cell(layer, m_farPositive);
      const Cell negativeEnd = m_geo.get_cell(layer, m_farNegative);
      return perSubPos(
          [&](Cell::SubPos subpos)
          {
            return Cylinder(std::abs(cell.r(subpos)),
                            0.5
                                * (std::abs(positiveEnd.z(subpos))
                                   + std::abs(negativeEnd.z(subpos))));
          });
    }

    // Out to the radius at the layer's smallest |eta| on the particle's side
    const double eta = m_positiveSide
        ? m_geo.min_eta(layer, CaloGeo::DetectorSide::kEtaPositive)
        : m_geo.max_eta(layer, CaloGeo::DetectorSide::kEtaNegative);
    const double theta = 2 * std::atan(std::exp(-eta));
    return perSubPos(
        [&](Cell::SubPos subpos)
        {
          const auto halfLength = static_cast<float>(std::abs(cell.z(subpos)));
          const double radius = std::abs(
              halfLength
              * std::sqrt((1 / (std::cos(theta) * std::cos(theta))) - 1));
          return Cylinder(radius, halfLength);
        });
  }

private:
  /// The z of the front of EME1 on the particle's side, looked up once
  auto embHalfLength() -> double
  {
    if (!m_embHalfLength) {
      m_embHalfLength =
          std::abs(m_geo.zpos(5,
                              m_positiveSide ? m_farPositive : m_farNegative,
                              Cell::SubPos::ENT));
    }
    return *m_embHalfLength;
  }

  const CaloGeo& m_geo;
  /// The position on the ID-calo boundary
  Position m_boundary;
  /// Positions far out on either side of the detector, at the particle's phi
  Position m_farPositive;
  Position m_farNegative;
  bool m_positiveSide;
  std::optional<double> m_embHalfLength;
};

/// Scales the position, keeping its eta and phi, onto the part of the cylinder
/// the layer is on: the cover for barrel layers, the endcap for the others
auto scaleOntoLayer(const Vector3D& position,
                    const Cylinder& cylinder,
                    bool isBarrel) -> Vector3D
{
  double scale = 1;
  if (isBarrel && position.perp() > 1e-6) {
    scale = cylinder.radius / position.perp();
  } else if (!isBarrel && std::abs(position.z()) > 1e-6) {
    scale = cylinder.halfLength / std::abs(position.z());
  }
  return scale * position;
}

auto describe(const TFCSTruthState& truth) -> std::string
{
  return fmt::format(
      "particle with PDG ID {}, pT={}, eta={}, phi={}, E={}, Ekin_off={} and "
      "vertex ({}, {}, {})",
      truth.pdgid(),
      truth.Pt(),
      truth.Eta(),
      truth.Phi(),
      truth.E(),
      truth.Ekin_off(),
      truth.vertex().X(),
      truth.vertex().Y(),
      truth.vertex().Z());
}

}  // namespace

void FastCaloSimCaloExtrapolation::extrapolate(
    TFCSExtrapolationState& result,
    const TFCSTruthState* truth,
    const std::vector<G4FieldTrack>& caloSteps) const
{
  const auto reportFailure = [&](const std::string& what)
  {
    const std::string message = what + " for " + describe(*truth);
    if (truth->Pt() < kWarningPtLimit) {
      FCS_MSG_DEBUG(message);
    } else {
      FCS_MSG_WARNING(message);
    }
  };

  std::vector<TrackPoint> track;
  track.reserve(caloSteps.size());
  for (const G4FieldTrack& step : caloSteps) {
    track.push_back({step.GetPosition(), step.GetMomentumDir()});
  }

  // ID-calo boundary
  const std::optional<TrackPoint> boundary = extrapolateToIDCaloBoundary(track);
  if (!boundary) {
    result.set_IDCaloBoundary_eta(-999.);
    result.set_IDCaloBoundary_phi(-999.);
    result.set_IDCaloBoundary_r(0);
    result.set_IDCaloBoundary_x(0);
    result.set_IDCaloBoundary_y(0);
    result.set_IDCaloBoundary_z(0);
    result.set_IDCaloBoundary_AngleEta(-999.);
    result.set_IDCaloBoundary_Angle3D(-999.);
    reportFailure("Failed extrapolation to the ID-calo boundary");
    return;
  }
  const Vector3D& position = boundary->position;
  const Vector3D& direction = boundary->direction;
  result.set_IDCaloBoundary_eta(position.eta());
  result.set_IDCaloBoundary_phi(position.phi());
  result.set_IDCaloBoundary_r(position.perp());
  result.set_IDCaloBoundary_x(position.x());
  result.set_IDCaloBoundary_y(position.y());
  result.set_IDCaloBoundary_z(position.z());
  // The angles between position and momentum, to correct the shower shapes of
  // particles that do not come from the origin
  result.set_IDCaloBoundary_AngleEta(position.theta() - direction.theta());
  result.set_IDCaloBoundary_Angle3D(position.angle(direction));
  FCS_MSG_DEBUG("ID-calo boundary at eta="
                << position.eta() << " phi=" << position.phi()
                << " r=" << position.perp() << " z=" << position.z());

  const bool reachesCalorimeter = std::abs(position.eta()) < 6;
  if (!reachesCalorimeter) {
    reportFailure(fmt::format(
        "ID-calo boundary at eta={} outside the calorimeter", position.eta()));
    return;
  }

  // Calorimeter layers
  LayerCylinders layerCylinders(*m_geo, position);
  for (unsigned int layer = 0; layer < m_geo->n_layers(); ++layer) {
    const bool isBarrel = m_geo->is_barrel(layer);
    const std::array<Cylinder, 3> cylinders = layerCylinders.of(layer);
    for (const auto subpos :
         {Cell::SubPos::MID, Cell::SubPos::ENT, Cell::SubPos::EXT})
    {
      const Cylinder& cylinder = cylinders[subpos];
      const auto point = extrapolateToCylinder(track, cylinder);
      if (!point) {
        reportFailure(fmt::format("Failed extrapolation to layer {} subpos {}",
                                  layer,
                                  static_cast<int>(subpos)));
        continue;
      }
      const Vector3D onLayer =
          scaleOntoLayer(point->position, cylinder, isBarrel);
      result.set_OK(layer, subpos, true);
      result.set_eta(layer, subpos, onLayer.eta());
      result.set_phi(layer, subpos, onLayer.phi());
      result.set_r(layer, subpos, onLayer.perp());
      result.set_z(layer, subpos, onLayer.z());
      FCS_MSG_DEBUG("Layer "
                    << layer << " subpos " << subpos
                    << " at eta=" << onLayer.eta() << " phi=" << onLayer.phi()
                    << " r=" << onLayer.perp() << " z=" << onLayer.z());
    }
  }
}
