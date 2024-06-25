#pragma once

#include <vector>

#include <nlohmann/json.hpp>

#include "CLHEP/Units/SystemOfUnits.h"
#include "TestHelpers/Container.h"
#include "Track.h"

namespace TestHelpers
{

struct TrackContainerData
{
  int track_id;
  double x, y, z, r, eta, phi;

  NLOHMANN_DEFINE_TYPE_INTRUSIVE(
      TrackContainerData, track_id, x, y, z, r, eta, phi)
};

class TrackContainer : public Container<Track, TrackContainerData>
{
  auto flatten() const -> std::vector<TrackContainerData> override
  {
    std::vector<TrackContainerData> data;
    int track_id = 0;
    for (const auto& track : m_container) {
      for (const auto& point : track) {
        data.push_back({track_id,
                             point.x() / CLHEP::m,
                             point.y() / CLHEP::m,
                             point.z() / CLHEP::m,
                             point.rho() / CLHEP::m,
                             point.eta(),
                             point.phi()});
      }
      ++track_id;
    }
    return data;
  }
};

}  // namespace TestHelpers
