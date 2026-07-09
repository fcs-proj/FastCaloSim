// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

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
  FCS_MSG_INFO("IDCalo: eta="
               << m_IDCaloBoundary_eta << " phi=" << m_IDCaloBoundary_phi
               << " r=" << m_IDCaloBoundary_r << " z=" << m_IDCaloBoundary_z);

  // Iterate over the unordered_map to print layer/subpos info
  for (const auto& [key, ok] : m_CaloOK) {
    int layer = key.first;
    int subpos = key.second;
    if (ok) {
      FCS_MSG_INFO("  layer " << layer << " subpos " << subpos
                              << " MID eta=" << m_etaCalo.at({layer, subpos})
                              << " phi=" << m_phiCalo.at({layer, subpos})
                              << " r=" << m_rCalo.at({layer, subpos})
                              << " z=" << m_zCalo.at({layer, subpos}));
    }
  }
}

void TFCSExtrapolationState::clear()
{
  // Drop all stored layer/subpos entries. The getters return the default
  // values (OK=false, eta/phi=-999, r/z/d/detaBorder=0) for entries that are
  // not stored, so no pre-filling is needed.
  m_CaloOK.clear();
  m_etaCalo.clear();
  m_phiCalo.clear();
  m_rCalo.clear();
  m_zCalo.clear();
  m_dCalo.clear();
  m_distetaCaloBorder.clear();

  // Reset IDCaloBoundary variables
  m_IDCaloBoundary_eta = -999;
  m_IDCaloBoundary_phi = -999;
  m_IDCaloBoundary_r = 0;
  m_IDCaloBoundary_z = 0;

  m_IDCaloBoundary_AngleEta = -999;
  m_IDCaloBoundary_Angle3D = -999;
}
