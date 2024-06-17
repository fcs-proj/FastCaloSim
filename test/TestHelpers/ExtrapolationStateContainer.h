#pragma once

#include <fstream>
#include <string>
#include <vector>

#include <CLHEP/Units/SystemOfUnits.h>
#include <nlohmann/json.hpp>

#include "FastCaloSim/Core/TFCSExtrapolationState.h"

namespace TestHelpers
{

struct ExtrapolationStateContainerData
{
  int state_id;
  int layer_id;
  std::string subpos;
  double eta, phi, r, z;

  NLOHMANN_DEFINE_TYPE_INTRUSIVE(ExtrapolationStateContainerData,
                                 state_id,
                                 layer_id,
                                 subpos,
                                 eta,
                                 phi,
                                 r,
                                 z)
};

class ExtrapolationStateContainer
{
public:
  ExtrapolationStateContainer() = default;

  explicit ExtrapolationStateContainer(
      std::vector<TFCSExtrapolationState> container)
      : m_container(std::move(container))
  {
  }

  /**
   * @brief Add an ExtrapolationState to the ExtrapolationState container
   */
  void addState(const TFCSExtrapolationState& state)
  {
    m_container.push_back(state);
  }

  /**
   * @brief Clear the ExtrapolationState container
   */
  void clear() { m_container.clear(); }

  /**
   * @brief Serialize the ExtrapolationState container to a JSON file
   * Distance units are converted to meters
   */
  void serialize(const std::string& filename)
  {
    std::ofstream ofs(filename);

    std::vector<ExtrapolationStateContainerData> positions;
    int state_id = 0;
    for (const auto& state : m_container) {
      // Serialize IDCaloBoundary with layer_id = -1
      positions.push_back({state_id,
                           -1,
                           "IDCaloBoundary",
                           state.IDCaloBoundary_eta(),
                           state.IDCaloBoundary_phi(),
                           state.IDCaloBoundary_r() / CLHEP::m,
                           state.IDCaloBoundary_z() / CLHEP::m});

      // Define subpositions
      const std::vector<std::pair<std::string, TFCSExtrapolationState::SUBPOS>>
          subpositions = {{"ENT", TFCSExtrapolationState::SUBPOS::SUBPOS_ENT},
                          {"MID", TFCSExtrapolationState::SUBPOS::SUBPOS_MID},
                          {"EXT", TFCSExtrapolationState::SUBPOS::SUBPOS_EXT}};

      // Serialize each layer for each subposition
      for (int layer = 0; layer <= 23; ++layer) {
        for (const auto& [subpos_name, subpos] : subpositions) {
          positions.push_back({state_id,
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

    nlohmann::json j = positions;
    ofs << j.dump();
    ofs.close();
  }

private:
  std::vector<TFCSExtrapolationState> m_container;
};

}  // namespace TestHelpers
