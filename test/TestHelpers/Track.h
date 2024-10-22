// Copyright (c) 2024 CERN for the benefit of the FastCaloSim project

#pragma once

#include <vector>

#include "CLHEP/Vector/ThreeVector.h"
#include "G4FieldTrack.hh"

namespace TestHelpers
{

class Track
{
  /**
   * @brief A track defined as a collection of points in 3D space
   */

  using Vector3D = CLHEP::Hep3Vector;

public:
  // Default constructor
  Track() = default;

  explicit Track(std::vector<Vector3D> positions)
      : m_track(std::move(positions))
  {
  }

  explicit Track(std::vector<G4FieldTrack> caloSteps)
  {
    for (const auto& step : caloSteps) {
      m_track.emplace_back(step.GetPosition());
    }
  }

  auto begin() const { return m_track.begin(); }
  auto end() const { return m_track.end(); }

private:
  std::vector<Vector3D> m_track;
};

}  // namespace TestHelpers
