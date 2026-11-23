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
  // Indexed by layer ID
  std::vector<double> energies(m_geo->layer_id_bound(), 0);

  // Loop over all cells containing energy and sum up energies
  for (const auto& cell_iter : simulstate.cells()) {
    unsigned long long const cell_id = cell_iter.first;
    const double cell_energy = cell_iter.second;
    // Add the energy to the layer of the cell
    const auto layer = m_geo->get_cell(cell_id).layer();
    energies[layer] += cell_energy;
  }

  std::vector<double> scalefactor(m_geo->layer_id_bound(), 1);

  // Expected layer energies above these thresholds that no cell received are
  // reported (none if the host set no layer roles)
  const std::vector<double>* report_threshold = m_geo->has_layer_roles()
      ? &m_geo->layer_roles().out_of_calo_report_threshold
      : nullptr;

  for (const std::size_t layer : m_geo->layer_ids()) {
    // catch large amounts of energy not simulated as shower is outside the
    // calorimeter
    if (FastCaloSim::exactly_equal(energies[layer], 0.0)
        && !FastCaloSim::exactly_equal(simulstate.E(layer), 0.0)
        && report_threshold != nullptr && layer < report_threshold->size()
        && simulstate.E(layer) > (*report_threshold)[layer])
    {
      FCS_MSG_INFO(
          "TFCSEnergyRenormalization::simulate(): energy not simulated "
          "(out-of-calo) in layer "
          << layer << " expected: " << simulstate.E(layer)
          << " simulated: " << energies[layer]);
    }
    if (!FastCaloSim::exactly_equal(energies[layer], 0.0)) {
      scalefactor[layer] = simulstate.E(layer) / energies[layer];
    }
  }

  // Apply the scale factor of its layer to each cell
  simulstate.scale_cells(
      [this, &scalefactor](unsigned long long cell_id) -> double
      { return scalefactor[m_geo->get_cell(cell_id).layer()]; });

  if (msgLvl(FastCaloSim::MSG::DEBUG)) {
    FCS_MSG_DEBUG("Apply scale factors : ");
    for (const std::size_t layer : m_geo->layer_ids()) {
      FCS_MSG_DEBUG("  " << layer << " *= " << scalefactor[layer] << " : "
                         << energies[layer] << " -> " << simulstate.E(layer));
    }
  }

  return FCSSuccess;
}
