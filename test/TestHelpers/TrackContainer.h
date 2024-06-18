#pragma once

#include <fstream>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

#include "CLHEP/Units/SystemOfUnits.h"
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

class TrackContainer
{
public:
  TrackContainer() = default;

  explicit TrackContainer(std::vector<TestHelpers::Track> container)
      : m_trackContainer(std::move(container))
  {
  }

  /**
   * @brief Add a track to the track container
   */
  void addTrack(const TestHelpers::Track& track)
  {
    m_trackContainer.push_back(track);
  }

  /**
   * @brief Clear the track container
   */
  void clear() { m_trackContainer.clear(); }

  /**
   * @brief Serialize the track container to a pickle file
   * Distance units are converted to meters
   */
  void serialize(const std::string& filename)
  {
    std::ofstream ofs(filename);

    std::vector<TrackContainerData> positions;
    int track_id = 0;
    for (const auto& track : m_trackContainer) {
      for (const auto& point : track) {
        positions.push_back({track_id,
                             point.x() / CLHEP::m,
                             point.y() / CLHEP::m,
                             point.z() / CLHEP::m,
                             point.rho() / CLHEP::m,
                             point.eta(),
                             point.phi()});
      }
      ++track_id;
    }

    nlohmann::json j = positions;
    ofs << j.dump();
    ofs.close();
  }

private:
  std::vector<TestHelpers::Track> m_trackContainer;
};

}  // namespace TestHelpers
