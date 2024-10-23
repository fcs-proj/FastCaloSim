// Copyright (c) 2024 CERN for the benefit of the FastCaloSim project

#include <iostream>

#include "FastCaloSim/Core/TFCSExtrapolationState.h"

//=============================================
//======= TFCSExtrapolationState =========
//=============================================

TFCSExtrapolationState::TFCSExtrapolationState()
{
  clear();
}

void TFCSExtrapolationState::Print(Option_t*) const
{
  // Print the IDCaloBoundary information
  ATH_MSG_INFO("IDCalo: eta="
               << m_IDCaloBoundary_eta << " phi=" << m_IDCaloBoundary_phi
               << " r=" << m_IDCaloBoundary_r << " z=" << m_IDCaloBoundary_z);

  // Iterate over the unordered_map to print layer/subpos info
  for (const auto& [key, ok] : m_CaloOK) {
    int layer = key.first;
    int subpos = key.second;
    if (ok) {
      ATH_MSG_INFO("  layer " << layer << " subpos " << subpos
                              << " MID eta=" << m_etaCalo.at({layer, subpos})
                              << " phi=" << m_phiCalo.at({layer, subpos})
                              << " r=" << m_rCalo.at({layer, subpos})
                              << " z=" << m_zCalo.at({layer, subpos}));
    }
  }
}

void TFCSExtrapolationState::clear()
{
  // Clear all unordered_maps and reset with default values
  for (int layer = 0; layer < m_CaloOK.size(); ++layer) {
    for (int subpos = 0; subpos < NumSubPos; ++subpos) {
      m_CaloOK[{layer, subpos}] = false;
      m_etaCalo[{layer, subpos}] = -999;
      m_phiCalo[{layer, subpos}] = -999;
      m_rCalo[{layer, subpos}] = 0;
      m_zCalo[{layer, subpos}] = 0;
      m_dCalo[{layer, subpos}] = 0;
      m_distetaCaloBorder[{layer, subpos}] = 0;
    }
  }

  // Reset IDCaloBoundary variables
  m_IDCaloBoundary_eta = -999;
  m_IDCaloBoundary_phi = -999;
  m_IDCaloBoundary_r = 0;
  m_IDCaloBoundary_z = 0;

  m_IDCaloBoundary_AngleEta = -999;
  m_IDCaloBoundary_Angle3D = -999;
}
