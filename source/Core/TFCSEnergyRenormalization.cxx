// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#include <cstddef>
#include <vector>

#include "FastCaloSim/Core/TFCSEnergyRenormalization.h"

#include "FastCaloSim/Core/TFCSSimulationState.h"
#include "FastCaloSim/Definitions/FloatCompare.h"
#include "FastCaloSim/Geometry/CaloGeo.h"
#include "FastCaloSim/Geometry/Cell.h"
#include "FastCaloSim/Geometry/LayerRoles.h"

//=============================================
//======= TFCSEnergyRenormalization =========
//=============================================

TFCSEnergyRenormalization::TFCSEnergyRenormalization(const char* name,
                                                     const char* title,
                                                     CaloGeo* geo)
    : TFCSParametrization(name, title)
    , m_geo(geo)
{
}

TFCSEnergyRenormalization::~TFCSEnergyRenormalization() = default;

auto TFCSEnergyRenormalization::simulate(
    TFCSSimulationState& simulstate,
    const TFCSTruthState* /*truth*/,
    const TFCSExtrapolationState* /*extrapol*/) const -> FCSReturnCode
{
  // One entry per layer, at its position in the layer IDs of the geometry
  const std::vector<std::size_t>& layers = m_geo->layer_ids();
  std::vector<double> energies(layers.size(), 0);

  // Loop over all cells containing energy and sum up energies
  for (const auto& cell_iter : simulstate.cells()) {
    unsigned long long const cell_id = cell_iter.first;
    const double cell_energy = cell_iter.second;
    // Add the energy to the layer of the cell
    energies[m_geo->layer_index(m_geo->get_cell(cell_id).layer())] +=
        cell_energy;
  }

  std::vector<double> scalefactor(layers.size(), 1);

  // Expected layer energies above these thresholds that no cell received are
  // reported (none if the host set no layer roles)
  const std::vector<double>* report_threshold = m_geo->has_layer_roles()
      ? &m_geo->layer_roles().out_of_calo_report_threshold
      : nullptr;

  for (std::size_t index = 0; index < layers.size(); ++index) {
    const std::size_t layer = layers[index];
    // catch large amounts of energy not simulated as shower is outside the
    // calorimeter
    if (FastCaloSim::exactly_equal(energies[index], 0.0)
        && !FastCaloSim::exactly_equal(simulstate.E(layer), 0.0)
        && report_threshold != nullptr && layer < report_threshold->size()
        && simulstate.E(layer) > (*report_threshold)[layer])
    {
      FCS_MSG_INFO(
          "TFCSEnergyRenormalization::simulate(): energy not simulated "
          "(out-of-calo) in layer "
          << layer << " expected: " << simulstate.E(layer)
          << " simulated: " << energies[index]);
    }
    if (!FastCaloSim::exactly_equal(energies[index], 0.0)) {
      scalefactor[index] = simulstate.E(layer) / energies[index];
    }
  }

  // Apply the scale factor of its layer to each cell
  simulstate.scale_cells(
      [this, &scalefactor](unsigned long long cell_id) -> double
      {
        return scalefactor[m_geo->layer_index(
            m_geo->get_cell(cell_id).layer())];
      });

  if (msgLvl(FastCaloSim::MSG::DEBUG)) {
    FCS_MSG_DEBUG("Apply scale factors : ");
    for (std::size_t index = 0; index < layers.size(); ++index) {
      FCS_MSG_DEBUG("  " << layers[index] << " *= " << scalefactor[index]
                         << " : " << energies[index] << " -> "
                         << simulstate.E(layers[index]));
    }
  }

  return FCSSuccess;
}
