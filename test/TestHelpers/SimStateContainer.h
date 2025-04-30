// Copyright (c) 2024 CERN for the benefit of the FastCaloSim project

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
  unsigned int state_id;  // The simulation state to which the cell belongs
  unsigned long long cell_id;  // The unique cell identifier
  float cell_energy;  // The deposited energy in that cell
  unsigned int layer;  // The layer ID to which the cell belons
  float eta, phi, r, x, y, z;  // The center positions of the cell
  float deta, dphi, dr, dx, dy, dz;  // The dimensions of the cell
  bool isBarrel;  // Whether the cell belongs to a barrel or endcap layer
  bool isXYZ, isEtaPhiR, isEtaPhiZ, isRPhiZ;  // The coordinate system in which
                                              // cell dimensions are described

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
                                 isXYZ,
                                 isEtaPhiR,
                                 isEtaPhiZ,
                                 isRPhiZ)
};

class SimStateContainer
    : public Container<TFCSSimulationState, SimStateContainerData>
{
  auto flatten() const -> std::vector<SimStateContainerData> override
  {
    std::vector<SimStateContainerData> data;

    int state_id = 0;
    for (const auto& state : m_container) {
      const auto& cellmap = state.cells();

      // Order by cell id
      std::map<unsigned long long, TFCSSimulationState::CellInfo>
          ordered_cellmap(cellmap.begin(), cellmap.end());

      for (const auto& [cell_id, info] : ordered_cellmap) {
        const Cell cell = m_geo->get_cell(cell_id);

        data.push_back({state_id,         cell.id(),        info.energy,
                        info.layer,       cell.eta(),       cell.phi(),
                        cell.r(),         cell.x(),         cell.y(),
                        cell.z(),         cell.deta(),      cell.dphi(),
                        cell.dr(),        cell.dx(),        cell.dy(),
                        cell.dz(),        cell.isBarrel(),  cell.isXYZ(),
                        cell.isEtaPhiR(), cell.isEtaPhiZ(), cell.isRPhiZ()});

        ++state_id;
      }
    }

    return data;
  };
};

}  // namespace TestHelpers
