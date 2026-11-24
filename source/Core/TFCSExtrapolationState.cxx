// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#include <iostream>

#include "FastCaloSim/Core/TFCSExtrapolationState.h"

//=============================================
//======= TFCSExtrapolationState =========
//=============================================

TFCSExtrapolationState::TFCSExtrapolationState() = default;

void TFCSExtrapolationState::Print(Option_t* /*option*/) const
{
  // Print the IDCaloBoundary information
  FCS_MSG_INFO("IDCalo: eta="
               << m_IDCaloBoundary_eta << " phi=" << m_IDCaloBoundary_phi
               << " r=" << m_IDCaloBoundary_r << " z=" << m_IDCaloBoundary_z);

  for (std::size_t i = 0; i < m_positions.size(); ++i) {
    const LayerPosition& pos = m_positions[i];
    if (pos.ok) {
      FCS_MSG_INFO("  layer " << i / NumSubPos << " subpos " << i % NumSubPos
                              << " eta=" << pos.eta << " phi=" << pos.phi
                              << " r=" << pos.r << " z=" << pos.z);
    }
  }
}

void TFCSExtrapolationState::clear()
{
  // Drop all stored positions: the getters return the defaults of a
  // sub-position that was not reached for positions that are not stored
  m_positions.clear();

  // Reset IDCaloBoundary variables
  m_IDCaloBoundary_eta = -999;
  m_IDCaloBoundary_phi = -999;
  m_IDCaloBoundary_r = 0;
  m_IDCaloBoundary_x = 0;
  m_IDCaloBoundary_y = 0;
  m_IDCaloBoundary_z = 0;

  m_IDCaloBoundary_AngleEta = -999;
  m_IDCaloBoundary_Angle3D = -999;
}
