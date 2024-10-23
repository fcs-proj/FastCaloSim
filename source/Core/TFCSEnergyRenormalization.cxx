#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-equal"
// Copyright (c) 2024 CERN for the benefit of the FastCaloSim project

#include "FastCaloSim/Core/TFCSEnergyRenormalization.h"

#include "FastCaloSim/Core/TFCSSimulationState.h"
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

TFCSEnergyRenormalization::~TFCSEnergyRenormalization() {}

FCSReturnCode TFCSEnergyRenormalization::simulate(
    TFCSSimulationState& simulstate,
    const TFCSTruthState* /*truth*/,
    const TFCSExtrapolationState* /*extrapol*/) const
{
  std::vector<double> energies(m_geo->n_layers(), 0);

  // Loop over all cells containing energy and sum up energies
  for (const auto& cell_iter : simulstate.cells()) {
    long long cell_id = cell_iter.first;
    double cell_energy = cell_iter.second;
    // Retrieve the cell from the geometry by its identifier
    Cell cell = m_geo->get_cell(cell_id);

    // Add the energy to the corresponding layer
    int layer = cell.layer();
    energies[layer] += cell_energy;
  }

  std::vector<float> scalefactor(m_geo->n_layers(), 1);

  const std::map<int, float> approxLayerNoise {
      {0, 150.0},  {1, 40.0},   {2, 80.0},   {3, 40.0},   {4, 150.0},
      {5, 40.0},   {6, 80.0},   {7, 50.0},   {8, 400.0},  {9, 400.0},
      {10, 400.0}, {11, 400.0}, {12, 300.0}, {13, 150.0}, {14, 40.0},
      {15, 150.0}, {16, 40.0},  {17, 400.0}, {18, 300.0}, {19, 150.0},
      {20, 40.0},  {21, 400.0}, {22, 400.0}, {23, 400.0}};

  for (int layer = 0; layer < m_geo->n_layers(); ++layer) {
    // catch large amounts of energy not simulated as shower is outside the
    // calorimeter
    if (energies[layer] == 0 && simulstate.E(layer) != 0) {
      if (simulstate.E(layer) > 8.0 * approxLayerNoise.at(layer) && layer != 5
          && layer != 6 && layer != 7)
        ATH_MSG_INFO(
            "TFCSEnergyRenormalization::simulate(): energy not simulated "
            "(out-of-calo) in layer "
            << layer << " expected: " << simulstate.E(layer)
            << " simulated: " << energies[layer]);
      if (simulstate.E(layer) > 1500.0
          && (layer == 5 || layer == 6 || layer == 7))
        ATH_MSG_INFO(
            "TFCSEnergyRenormalization::simulate(): energy not simulated "
            "(out-of-calo) in layer "
            << layer << " expected: " << simulstate.E(layer)
            << " simulated: " << energies[layer]);
    }
    if (energies[layer] != 0)
      scalefactor[layer] = simulstate.E(layer) / energies[layer];
  }

  // Loop over all cells and apply the scalefactor
  for (auto& cell_iter : simulstate.cells()) {
    long long cell_id = cell_iter.first;
    double cell_energy = cell_iter.second;
    Cell cell = m_geo->get_cell(cell_id);
    int layer = cell.layer();
    cell_iter.second *= scalefactor[layer];
  }

  if (msgLvl(MSG::DEBUG)) {
    ATH_MSG_DEBUG("Apply scale factors : ");
    for (int layer = 0; layer < m_geo->n_layers(); ++layer) {
      ATH_MSG_DEBUG("  " << layer << " *= " << scalefactor[layer] << " : "
                         << energies[layer] << " -> " << simulstate.E(layer));
    }
  }

  return FCSSuccess;
}
#pragma GCC diagnostic pop
