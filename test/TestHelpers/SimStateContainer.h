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
  int state_id;  // The simulation state to which the cell bellongs
  unsigned long long cell_id;  // The unique cell identifier
  float cell_energy;  // The deposited energy in that cell
  int layer_id;  // The layer ID to which the cell belons
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
                                 layer_id,
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

        // Only support hard-coded layer 1 for now
        if (cell->getSampling() != 1)
          continue;

        data.push_back({
            state_id,
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
            cell->dr(),
            cell->dx(),
            cell->dy(),
            cell->dz(),
            // TODO: the following is currently hardcoded for layer 1
            // should be adapted once we load the geometry with the exp-ind file
            1,  // cell->isBarrel(),
            1,  // cell->isCylindrical(),
            0,  // cell->isECCylindrical(),
            0  // cell->isCartesian()});
        });

        ++state_id;
      }
    }

    return data;
  };
};

}  // namespace TestHelpers
