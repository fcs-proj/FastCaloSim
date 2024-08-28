#pragma once

#include <vector>

#include <CLHEP/Units/SystemOfUnits.h>
#include <nlohmann/json.hpp>

#include "FastCaloSim/Core/TFCSSimulationState.h"
#include "FastCaloSim/Geometry/CaloDetDescrElement.h"
#include "TestHelpers/Container.h"

namespace TestHelpers
{

struct SimStateContainerData
{
  int state_id;  // The simulation state to which the cell belongs
  unsigned long long cell_id;  // The unique cell identifier
  float cell_energy;  // The deposited energy in that cell
  int layer;  // The layer ID to which the cell belons
  float eta, phi, r, x, y, z;  // The center positions of the cell
  float deta, dphi, dr, dx, dy, dz;  // The dimensions of the cell
  bool isBarrel;  // Whether the cell belongs to a barrel or endcap layer
  bool isCylindrical, isECCylindrical,
      isCartesian;  // The coordinate system in which cell dimensions are
                    // described

  NLOHMANN_DEFINE_TYPE_INTRUSIVE(SimStateContainerData,
                                 state_id,
                                 cell_id,
                                 cell_energy,
                                 layer,
                                 eta,
                                 phi,
                                 r,
                                 x,
                                 y,
                                 z,
                                 deta,
                                 dphi,
                                 dr,
                                 dx,
                                 dy,
                                 dz,
                                 isBarrel,
                                 isCylindrical,
                                 isECCylindrical,
                                 isCartesian)
};

class SimStateContainer
    : public Container<TFCSSimulationState, SimStateContainerData>
{
  auto flatten() const -> std::vector<SimStateContainerData> override
  {
    std::vector<SimStateContainerData> data;

    int state_id = 0;
    for (const auto& state : m_container) {
      TFCSSimulationState::Cellmap_t cellmap = state.cells();
      for (const auto& cell_map_element : cellmap) {
        const CaloDetDescrElement* cell = cell_map_element.first;
        const float cell_energy = cell_map_element.second;

        data.push_back(
            {state_id,
             cell->calo_hash(),
             cell_energy,
             cell->getSampling(),
             cell->eta(),
             cell->phi(),
             cell->r(),
             cell->x(),
             cell->y(),
             cell->z(),
             cell->deta(),
             cell->dphi(),
             cell->getSampling() < 4
                 ? cell->dr() * 2
                 : cell->dr(),  // dr saved here is only half cell size for
                                // barrel layers 0-3
             cell->dx(),
             cell->dy(),
             cell->dz() * 2,  // dz saved here is only half cell size
             // TODO: the following is currently hardcoded for ATLAS
             // should be adapted once we load the geometry with the exp-ind
             // file
             is_layer_type(cell->getSampling(), LayerType::BARREL),
             is_layer_type(cell->getSampling(), LayerType::CYLINDRICAL),
             is_layer_type(cell->getSampling(), LayerType::ECCYLINDRICAL),
             is_layer_type(cell->getSampling(), LayerType::CARTESIAN)});

        ++state_id;
      }
    }

    return data;
  };

  // hardcoded for now, this can be removed and adapted once we load the
  // geometry with the exp-ind file
  enum LayerType
  {
    BARREL,
    CYLINDRICAL,
    ECCYLINDRICAL,
    CARTESIAN
  };

  static auto is_layer_type(int sampling, LayerType type) -> bool
  {
    static const std::unordered_map<LayerType, std::vector<int>> layer_map = {
        {BARREL, {0, 1, 2, 3, 12, 13, 14, 15, 16, 17, 18, 19, 20}},
        {CYLINDRICAL, {0, 1, 2, 3, 12, 13, 14, 15, 16, 17, 18, 19, 20}},
        {ECCYLINDRICAL, {4, 5, 6, 7, 8, 9, 10, 11}},
        {CARTESIAN, {21, 22, 23}}};

    auto it = layer_map.find(type);
    if (it != layer_map.end()) {
      const std::vector<int>& layers = it->second;
      return std::find(layers.begin(), layers.end(), sampling) != layers.end();
    }
    return false;
  }
};

}  // namespace TestHelpers
