// Copyright (c) 2024 CERN for the benefit of the FastCaloSim project

#pragma once

#include <string>
#include <vector>

#include <CLHEP/Units/SystemOfUnits.h>
#include <nlohmann/json.hpp>

#include "FastCaloSim/Core/TFCSExtrapolationState.h"
#include "TestHelpers/Container.h"

namespace TestHelpers
{

struct ExtrapStateContainerData
{
  int state_id;
  int layer_id;
  std::string subpos;
  double eta, phi, r, z;

  NLOHMANN_DEFINE_TYPE_INTRUSIVE(
      ExtrapStateContainerData, state_id, layer_id, subpos, eta, phi, r, z)
};

class ExtrapStateContainer
    : public Container<TFCSExtrapolationState, ExtrapStateContainerData>
{
  auto flatten() const -> std::vector<ExtrapStateContainerData> override
  {
    std::vector<ExtrapStateContainerData> data;
    int state_id = 0;
    for (const auto& state : m_container) {
      // Serialize IDCaloBoundary with layer_id = -1
      data.push_back({state_id,
                      -1,
                      "IDCaloBoundary",
                      state.IDCaloBoundary_eta(),
                      state.IDCaloBoundary_phi(),
                      state.IDCaloBoundary_r() / CLHEP::m,
                      state.IDCaloBoundary_z() / CLHEP::m});

      // Define subpositions
      const std::vector<std::pair<std::string, Cell::SubPos>> subpositions = {
          {"ENT", Cell::SubPos::ENT},
          {"MID", Cell::SubPos::MID},
          {"EXT", Cell::SubPos::EXT}};

      // Serialize each layer for each subposition
      for (int layer = 0; layer <= 23; ++layer) {
        for (const auto& [subpos_name, subpos] : subpositions) {
          data.push_back({state_id,
                          layer,
                          subpos_name,
                          state.eta(layer, subpos),
                          state.phi(layer, subpos),
                          state.r(layer, subpos) / CLHEP::m,
                          state.z(layer, subpos) / CLHEP::m});
        }
      }

      ++state_id;
    }

    return data;
  }
};

}  // namespace TestHelpers
