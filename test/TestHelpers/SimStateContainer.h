#pragma once

#include <vector>

#include <CLHEP/Units/SystemOfUnits.h>
#include <nlohmann/json.hpp>

#include "FastCaloSim/Core/TFCSSimulationState.h"
#include "FastCaloSim/Geometry/CaloGeo.h"
#include "TestHelpers/Container.h"

namespace TestHelpers
{

struct SimStateContainerData
{
  int state_id;  // The simulation state to which the cell belongs
  long long cell_id;  // The unique cell identifier
  double cell_energy;  // The deposited energy in that cell
  long long layer;  // The layer ID to which the cell belons
  double eta, phi, r, x, y, z;  // The center positions of the cell
  double deta, dphi, dr, dx, dy, dz;  // The dimensions of the cell
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
      for (const auto& cell_map_iter : cellmap) {
        long long cell_id = cell_map_iter.first;
        const float cell_energy = cell_map_iter.second;

        // Get the cell
        const Cell cell = m_geo->get_cell(cell_id);

        // Get the cell dimensions where they make sense
        double dx = cell.isXYZ() ? cell.dx() : 0;
        double dy = cell.isXYZ() ? cell.dy() : 0;
        double dz = cell.isXYZ() || cell.isEtaPhiZ() ? cell.dz() : 0;
        double deta = cell.isEtaPhiR() || cell.isEtaPhiZ() ? cell.deta() : 0;
        double dphi = cell.isEtaPhiR() || cell.isEtaPhiZ() ? cell.dphi() : 0;
        double dr = cell.isEtaPhiR() ? cell.dr() : 0;

        data.push_back({state_id,
                        cell.id(),
                        cell_energy,
                        cell.layer(),
                        cell.eta(),
                        cell.phi(),
                        cell.r(),
                        cell.x(),
                        cell.y(),
                        cell.z(),
                        deta,
                        dphi,
                        dr,
                        dx,
                        dy,
                        dz,
                        cell.isBarrel(),
                        cell.isEtaPhiR(),
                        cell.isEtaPhiZ(),
                        cell.isXYZ()});

        ++state_id;
      }
    }

    return data;
  };
};

}  // namespace TestHelpers
