#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-equal"
/*
  Copyright (C) 2002-2019 CERN for the benefit of the ATLAS collaboration
*/

#include "FastCaloSim/Core/TFCSEnergyRenormalization.h"

#include "FastCaloSim/Geometry/CaloDetDescrElement.h"
#include "FastCaloSim/Geometry/FastCaloSim_CaloCell_ID.h"
#include "FastCaloSim/Core/TFCSSimulationState.h"

//=============================================
//======= TFCSEnergyRenormalization =========
//=============================================

TFCSEnergyRenormalization::TFCSEnergyRenormalization(const char* name,
                                                     const char* title)
    : TFCSParametrization(name, title)
{
}

TFCSEnergyRenormalization::~TFCSEnergyRenormalization() {}

FCSReturnCode TFCSEnergyRenormalization::simulate(
    TFCSSimulationState& simulstate,
    const TFCSTruthState* /*truth*/,
    const TFCSExtrapolationState* /*extrapol*/) const
{
  std::vector<double> energies(CaloCell_ID_FCS::MaxSample, 0);

  // Loop over all cells and sum up energies
  for (const auto& iter : simulstate.cells()) {
    const CaloDetDescrElement* theDDE = iter.first;
    int layer = theDDE->getSampling();
    energies[layer] += iter.second;
  }

  std::vector<float> scalefactor(CaloCell_ID_FCS::MaxSample, 1);

  const std::map<int, float> approxLayerNoise {
      {0, 150.0},  {1, 40.0},   {2, 80.0},   {3, 40.0},   {4, 150.0},
      {5, 40.0},   {6, 80.0},   {7, 50.0},   {8, 400.0},  {9, 400.0},
      {10, 400.0}, {11, 400.0}, {12, 300.0}, {13, 150.0}, {14, 40.0},
      {15, 150.0}, {16, 40.0},  {17, 400.0}, {18, 300.0}, {19, 150.0},
      {20, 40.0},  {21, 400.0}, {22, 400.0}, {23, 400.0}};

  for (int layer = 0; layer < CaloCell_ID_FCS::MaxSample; ++layer) {
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
  for (auto& iter : simulstate.cells()) {
    const CaloDetDescrElement* theDDE = iter.first;
    int layer = theDDE->getSampling();
    iter.second *= scalefactor[layer];
  }

  if (msgLvl(MSG::DEBUG)) {
    ATH_MSG_DEBUG("Apply scale factors : ");
    for (int layer = 0; layer < CaloCell_ID_FCS::MaxSample; ++layer) {
      ATH_MSG_DEBUG("  " << layer << " *= " << scalefactor[layer] << " : "
                         << energies[layer] << " -> " << simulstate.E(layer));
    }
  }

  return FCSSuccess;
}
#pragma GCC diagnostic pop
