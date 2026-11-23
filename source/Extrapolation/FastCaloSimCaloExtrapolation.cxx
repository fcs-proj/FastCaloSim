// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <format>
#include <map>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "FastCaloSim/Extrapolation/FastCaloSimCaloExtrapolation.h"

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

using BoundaryCylinder = FastCaloSimCaloExtrapolation::BoundaryCylinder;

auto toCylinder(const BoundaryCylinder& cylinder) -> Cylinder
{
  return {cylinder.radius, cylinder.half_length};
}

/// Where the track leaves the ID: of the points where it reaches the boundary
/// cylinders, the closest to the origin that is on the boundary, i.e. neither
/// within the z range of the previous cylinder nor the radius of the next one
auto extrapolateToIDCaloBoundary(const std::vector<TrackPoint>& track,
                                 const std::vector<BoundaryCylinder>& cylinders)
    -> std::optional<TrackPoint>
{
  constexpr double tolerance = 0.001;  // mm

  std::optional<TrackPoint> boundary;
  for (std::size_t i = 0; i < cylinders.size(); ++i) {
    const auto point = extrapolateToCylinder(track, toCylinder(cylinders[i]));
    if (!point) {
      continue;
    }
    const Vector3D& position = point->position;
    // Up to the end of the previous cylinder, including the edge where the
    // two meet, which the previous cylinder's endcap already covers. A
    // cylinder that ends at the same |z| as the previous one closes the
    // previous one's endcap towards the beam axis: its own endcap is on the
    // boundary (a particle in the beam pipe reaches it there).
    const bool withinPrevious = i > 0
        && std::abs(position.z())
            < std::min(toCylinder(cylinders[i - 1]).half_length() + tolerance,
                       toCylinder(cylinders[i]).half_length() - tolerance);
    const bool withinNext = i + 1 < cylinders.size()
        && position.perp() < toCylinder(cylinders[i + 1]).r() + tolerance;
    if (withinPrevious || withinNext) {
      continue;
    }
    if (!boundary || position.mag() <= boundary->position.mag()) {
      boundary = point;
    }
  }
  return boundary;
}

/// The three cylinders of a layer, indexed by FastCaloSim::SubPos
template<typename CylinderOfSubPos>
auto perSubPos(const CylinderOfSubPos& cylinder) -> std::array<Cylinder, 3>
{
  static_assert(std::to_underlying(FastCaloSim::SubPos::MID) == 0
                && std::to_underlying(FastCaloSim::SubPos::ENT) == 1
                && std::to_underlying(FastCaloSim::SubPos::EXT) == 2);
  return {cylinder(FastCaloSim::SubPos::MID),
          cylinder(FastCaloSim::SubPos::ENT),
          cylinder(FastCaloSim::SubPos::EXT)};
}

/// The cylinders through the middle, entrance and exit of the calorimeter
/// layers, for a particle that crossed the ID-calo boundary at a position.
/// They go through the cell at that position: barrel layers are cylinder
/// covers, closed at the end of the layer or at the front of the layer given in
/// the settings; the other layers are endcaps, out to the edge of the layer
/// closest to the barrel.
class LayerCylinders
{
public:
  LayerCylinders(const CaloGeo& geo,
                 const std::map<std::size_t, std::size_t>& barrelLayerEnds,
                 const Vector3D& boundary)
      : m_geo(&geo)
      , m_barrelLayerEnds(&barrelLayerEnds)
      , m_boundary {.m_x = boundary.x(),
                    .m_y = boundary.y(),
                    .m_z = boundary.z(),
                    .m_eta = boundary.eta(),
                    .m_phi = boundary.phi(),
                    .m_r = boundary.perp()}
      , m_farPositive {.m_x = 0,
                       .m_y = 0,
                       .m_z = 0,
                       .m_eta = 1000,
                       .m_phi = boundary.phi(),
                       .m_r = 0}
      , m_farNegative {.m_x = 0,
                       .m_y = 0,
                       .m_z = 0,
                       .m_eta = -1000,
                       .m_phi = boundary.phi(),
                       .m_r = 0}
      , m_positiveSide(boundary.eta() > 0)
  {
  }

  /// The cylinders of the layer, indexed by FastCaloSim::SubPos
  auto of(std::size_t layer) -> std::array<Cylinder, 3>
  {
    const Cell cell = m_geo->get_cell(layer, m_boundary);

    if (m_geo->is_barrel(layer)) {
      if (const auto end = m_barrelLayerEnds->find(layer);
          end != m_barrelLayerEnds->end())
      {
        // Closed at the front of the other layer on the particle's side
        const double halfLength = frontOf(end->second);
        return perSubPos([&](FastCaloSim::SubPos subpos) -> Cylinder
                         { return {std::abs(cell.r(subpos)), halfLength}; });
      }
      // Closed at the mean |z| of the layer's ends on both sides
      const Cell positiveEnd = m_geo->get_cell(layer, m_farPositive);
      const Cell negativeEnd = m_geo->get_cell(layer, m_farNegative);
      return perSubPos(
          [&](FastCaloSim::SubPos subpos) -> Cylinder
          {
            return {std::abs(cell.r(subpos)),
                    0.5
                        * (std::abs(positiveEnd.z(subpos))
                           + std::abs(negativeEnd.z(subpos)))};
          });
    }

    // Out to the radius at the layer's smallest |eta| on the particle's side
    const double eta = m_positiveSide
        ? m_geo->min_eta(layer, CaloGeo::DetectorSide::kEtaPositive)
        : m_geo->max_eta(layer, CaloGeo::DetectorSide::kEtaNegative);
    const double theta = 2 * std::atan(std::exp(-eta));
    return perSubPos(
        [&](FastCaloSim::SubPos subpos) -> Cylinder
        {
          // The half length is rounded to single precision first
          const auto halfLength =
              static_cast<double>(static_cast<float>(std::abs(cell.z(subpos))));
          const double radius = std::abs(
              halfLength
              * std::sqrt((1 / (std::cos(theta) * std::cos(theta))) - 1));
          return {radius, halfLength};
        });
  }

private:
  /// The |z| of the front of a layer on the particle's side, looked up once
  auto frontOf(std::size_t layer) -> double
  {
    for (const auto& [frontLayer, halfLength] : m_fronts) {
      if (frontLayer == layer) {
        return halfLength;
      }
    }
    const double halfLength =
        std::abs(m_geo->zpos(layer,
                             m_positiveSide ? m_farPositive : m_farNegative,
                             FastCaloSim::SubPos::ENT));
    m_fronts.emplace_back(layer, halfLength);
    return halfLength;
  }

  const CaloGeo* m_geo;
  const std::map<std::size_t, std::size_t>* m_barrelLayerEnds;
  /// The position on the ID-calo boundary
  Position m_boundary;
  /// Positions far out on either side of the detector, at the particle's phi
  Position m_farPositive;
  Position m_farNegative;
  bool m_positiveSide;
  /// The fronts looked up so far: layer and |z|
  std::vector<std::pair<std::size_t, double>> m_fronts;
};

/// Scales the position, keeping its eta and phi, onto the part of the cylinder
/// the layer is on: the cover for barrel layers, the endcap for the others
auto scaleOntoLayer(const Vector3D& position,
                    const Cylinder& cylinder,
                    bool isBarrel) -> Vector3D
{
  double scale = 1;
  if (isBarrel && position.perp() > 1e-6) {
    scale = cylinder.r() / position.perp();
  } else if (!isBarrel && std::abs(position.z()) > 1e-6) {
    scale = cylinder.half_length() / std::abs(position.z());
  }
  return scale * position;
}

auto describe(const TFCSTruthState& truth) -> std::string
{
  return std::format(
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

/// Throws std::invalid_argument if the barrel layer ends of the settings do not
/// fit the geometry
void checkLayers(const CaloGeo& geo,
                 const FastCaloSimCaloExtrapolation::Settings& settings)
{
  const std::size_t n_layers = geo.n_layers();
  for (const auto& [barrelLayer, frontLayer] :
       settings.barrel_layer_ends_at_layer)
  {
    if (barrelLayer >= n_layers || !geo.is_barrel(barrelLayer)) {
      throw std::invalid_argument(std::format(
          "FastCaloSimCaloExtrapolation: layer {} is set to end at another "
          "layer, but is not a barrel layer of the geometry",
          barrelLayer));
    }
    if (frontLayer >= n_layers) {
      throw std::invalid_argument(std::format(
          "FastCaloSimCaloExtrapolation: barrel layer {} is set to end at "
          "layer {}, which the geometry does not have",
          barrelLayer,
          frontLayer));
    }
  }
}

}  // namespace

void FastCaloSimCaloExtrapolation::set_settings(Settings settings)
{
  const auto& cylinders = settings.id_calo_boundary;
  if (cylinders.empty()) {
    throw std::invalid_argument(
        "FastCaloSimCaloExtrapolation: the ID-calo boundary has no cylinders");
  }
  for (std::size_t i = 0; i < cylinders.size(); ++i) {
    const BoundaryCylinder& cylinder = cylinders[i];
    // Written so that NaN dimensions are rejected as well
    const bool positive = cylinder.radius > 0 && cylinder.half_length > 0;
    if (!positive) {
      throw std::invalid_argument(std::format(
          "FastCaloSimCaloExtrapolation: ID-calo boundary cylinder {} has "
          "radius {} and half length {}, both must be positive",
          i,
          cylinder.radius,
          cylinder.half_length));
    }
    const bool nested = i == 0
        || (cylinder.radius < cylinders[i - 1].radius
            && cylinder.half_length >= cylinders[i - 1].half_length);
    if (!nested) {
      throw std::invalid_argument(std::format(
          "FastCaloSimCaloExtrapolation: ID-calo boundary cylinder {} must be "
          "closer to the beam axis than the previous one and reach at least "
          "as far in |z|",
          i));
    }
  }
  if (m_geo != nullptr) {
    checkLayers(*m_geo, settings);
  }
  m_settings = std::move(settings);
}

void FastCaloSimCaloExtrapolation::set_geometry(CaloGeo* geo)
{
  if (geo != nullptr && !m_settings.id_calo_boundary.empty()) {
    checkLayers(*geo, m_settings);
  }
  m_geo = geo;
}

void FastCaloSimCaloExtrapolation::extrapolate(
    TFCSExtrapolationState& result,
    const TFCSTruthState* truth,
    const std::vector<G4FieldTrack>& caloSteps) const
{
  if (m_settings.id_calo_boundary.empty() || m_geo == nullptr) {
    throw std::logic_error(
        "FastCaloSimCaloExtrapolation: set_settings and set_geometry must be "
        "called before the first extrapolation");
  }

  const auto reportFailure = [&](const std::string& what) -> void
  {
    const std::string message = what + " for " + describe(*truth);
    if (truth->Pt() < m_settings.min_warning_pt) {
      FCS_MSG_DEBUG(message);
    } else {
      FCS_MSG_WARNING(message);
    }
  };

  std::vector<TrackPoint> track;
  track.reserve(caloSteps.size());
  for (const G4FieldTrack& step : caloSteps) {
    track.push_back(
        {.position = step.GetPosition(), .direction = step.GetMomentumDir()});
  }

  // ID-calo boundary
  const std::optional<TrackPoint> boundary =
      extrapolateToIDCaloBoundary(track, m_settings.id_calo_boundary);
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

  const bool reachesCalorimeter =
      std::abs(position.eta()) < m_settings.max_abs_eta;
  if (!reachesCalorimeter) {
    reportFailure(std::format(
        "ID-calo boundary at eta={} outside the calorimeter", position.eta()));
    return;
  }

  // Calorimeter layers
  LayerCylinders layerCylinders(
      *m_geo, m_settings.barrel_layer_ends_at_layer, position);
  for (std::size_t layer = 0; layer < m_geo->n_layers(); ++layer) {
    const bool isBarrel = m_geo->is_barrel(layer);
    const std::array<Cylinder, 3> cylinders = layerCylinders.of(layer);
    for (const auto subpos : {FastCaloSim::SubPos::MID,
                              FastCaloSim::SubPos::ENT,
                              FastCaloSim::SubPos::EXT})
    {
      const Cylinder& cylinder = cylinders[std::to_underlying(subpos)];
      const auto point = extrapolateToCylinder(track, cylinder);
      if (!point) {
        reportFailure(std::format("Failed extrapolation to layer {} subpos {}",
                                  layer,
                                  static_cast<int>(subpos)));
        continue;
      }
      const Vector3D onLayer =
          scaleOntoLayer(point->position, cylinder, isBarrel);
      result.set_OK(layer, subpos, /*val=*/true);
      result.set_eta(layer, subpos, onLayer.eta());
      result.set_phi(layer, subpos, onLayer.phi());
      result.set_r(layer, subpos, onLayer.perp());
      result.set_z(layer, subpos, onLayer.z());
      FCS_MSG_DEBUG("Layer "
                    << layer << " subpos " << static_cast<int>(subpos)
                    << " at eta=" << onLayer.eta() << " phi=" << onLayer.phi()
                    << " r=" << onLayer.perp() << " z=" << onLayer.z());
    }
  }
}
