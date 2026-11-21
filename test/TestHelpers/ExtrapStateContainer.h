// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

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
protected:
  auto flatten() const -> std::vector<ExtrapStateContainerData> override
  {
    std::vector<ExtrapStateContainerData> data;
    int state_id = 0;
    for (const auto& state : m_container) {
      // Serialize IDCaloBoundary with layer_id = -1
      data.push_back({.state_id = state_id,
                      .layer_id = -1,
                      .subpos = "IDCaloBoundary",
                      .eta = state.IDCaloBoundary_eta(),
                      .phi = state.IDCaloBoundary_phi(),
                      .r = state.IDCaloBoundary_r() / CLHEP::m,
                      .z = state.IDCaloBoundary_z() / CLHEP::m});

      // Define subpositions
      const std::vector<std::pair<std::string, FastCaloSim::SubPos>>
          subpositions = {{"ENT", FastCaloSim::SubPos::ENT},
                          {"MID", FastCaloSim::SubPos::MID},
                          {"EXT", FastCaloSim::SubPos::EXT}};

      // Serialize each layer for each subposition
      for (std::size_t layer = 0; layer <= 23; ++layer) {
        for (const auto& [subpos_name, subpos] : subpositions) {
          data.push_back({.state_id = state_id,
                          // The output column stays int
                          .layer_id = static_cast<int>(layer),
                          .subpos = subpos_name,
                          .eta = state.eta(layer, subpos),
                          .phi = state.phi(layer, subpos),
                          .r = state.r(layer, subpos) / CLHEP::m,
                          .z = state.z(layer, subpos) / CLHEP::m});
        }
      }

      ++state_id;
    }

    return data;
  }
};

}  // namespace TestHelpers
