// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#include <array>

#include "FastCaloSim/Core/TFCSEnergyRenormalization.h"

#include "FastCaloSim/Core/TFCSSimulationState.h"
#include "FastCaloSim/Definitions/FloatCompare.h"
#include "FastCaloSim/Geometry/CaloGeo.h"
#include "FastCaloSim/Geometry/Cell.h"

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
  std::vector<double> energies(m_geo->n_layers(), 0);

  // Loop over all cells containing energy and sum up energies
  for (const auto& cell_iter : simulstate.cells()) {
    unsigned long long const cell_id = cell_iter.first;
    const double cell_energy = cell_iter.second;
    // Add the energy to the layer of the cell
    const auto layer = m_geo->get_cell(cell_id).layer();
    energies[layer] += cell_energy;
  }

  std::vector<double> scalefactor(m_geo->n_layers(), 1);

  // Approximate noise per calorimeter layer in MeV
  static constexpr std::array<float, 24> approxLayerNoise {
      150.0, 40.0,  80.0,  40.0,  150.0, 40.0,  80.0,  50.0,
      400.0, 400.0, 400.0, 400.0, 300.0, 150.0, 40.0,  150.0,
      40.0,  400.0, 300.0, 150.0, 40.0,  400.0, 400.0, 400.0};

  for (std::size_t layer = 0; layer < m_geo->n_layers(); ++layer) {
    // catch large amounts of energy not simulated as shower is outside the
    // calorimeter
    if (FastCaloSim::exactly_equal(energies[layer], 0.0)
        && !FastCaloSim::exactly_equal(simulstate.E(layer), 0.0))
    {
      if (simulstate.E(layer)
              > 8.0 * static_cast<double>(approxLayerNoise.at(layer))
          && layer != 5 && layer != 6 && layer != 7)
      {
        FCS_MSG_INFO(
            "TFCSEnergyRenormalization::simulate(): energy not simulated "
            "(out-of-calo) in layer "
            << layer << " expected: " << simulstate.E(layer)
            << " simulated: " << energies[layer]);
      }
      if (simulstate.E(layer) > 1500.0
          && (layer == 5 || layer == 6 || layer == 7))
      {
        FCS_MSG_INFO(
            "TFCSEnergyRenormalization::simulate(): energy not simulated "
            "(out-of-calo) in layer "
            << layer << " expected: " << simulstate.E(layer)
            << " simulated: " << energies[layer]);
      }
    }
    if (!FastCaloSim::exactly_equal(energies[layer], 0.0)) {
      scalefactor[layer] = simulstate.E(layer) / energies[layer];
    }
  }

  // Loop over all cells and apply the scalefactor
  for (auto& cell_iter : simulstate.cells()) {
    const auto layer = m_geo->get_cell(cell_iter.first).layer();
    cell_iter.second *= scalefactor[layer];
  }

  if (msgLvl(FastCaloSim::MSG::DEBUG)) {
    FCS_MSG_DEBUG("Apply scale factors : ");
    for (std::size_t layer = 0; layer < m_geo->n_layers(); ++layer) {
      FCS_MSG_DEBUG("  " << layer << " *= " << scalefactor[layer] << " : "
                         << energies[layer] << " -> " << simulstate.E(layer));
    }
  }

  return FCSSuccess;
}
