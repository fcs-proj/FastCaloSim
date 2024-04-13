#pragma GCC diagnostic push 
#pragma GCC diagnostic ignored "-Wfloat-equal"

/*
  Copyright (C) 2002-2024 CERN for the benefit of the ATLAS collaboration
*/

#include "CLHEP/Random/RandomEngine.h"

#include "FastCaloSim/TFCSSimulationState.h"
#include "FastCaloSim/TFCSParametrizationBase.h"
#include <iostream>
#include <cstring>

//=============================================
//======= TFCSSimulationState =========
//=============================================

TFCSSimulationState::TFCSSimulationState(CLHEP::HepRandomEngine *randomEngine)
    : m_randomEngine(randomEngine) {
  clear();
}

void TFCSSimulationState::clear() {
  set_SF(1);
  m_Ebin = -1;
  m_Etot = 0;
  for (int i = 0; i < CaloCell_ID_FCS::MaxSample; ++i) {
    m_E[i] = 0;
    m_Efrac[i] = 0;
  }
}

void TFCSSimulationState::deposit(const CaloDetDescrElement *cellele, float E) {
  auto mele = m_cells.find(cellele);
  if (mele == m_cells.end()) {
    m_cells.emplace(cellele, 0);
    mele = m_cells.find(cellele);
  }
  m_cells[cellele] += E;
}

void TFCSSimulationState::Print(Option_t *) const {
  ATH_MSG_INFO("Ebin=" << m_Ebin << " E=" << E()
                       << " #cells=" << m_cells.size());
  for (int i = 0; i < CaloCell_ID_FCS::MaxSample; ++i)
    if (E(i) != 0) {
      ATH_MSG_INFO("  E" << i << "(" << CaloSampling::getSamplingName(i)
                         << ")=" << E(i) << " E" << i << "/E=" << Efrac(i));
    }
  if (!m_AuxInfo.empty()) {
    ATH_MSG_INFO("  AuxInfo has " << m_AuxInfo.size() << " elements");
    for (const auto &a : m_AuxInfo) {
      ATH_MSG_INFO("    " << a.first << " : "
                          // Dont print as char/bool.
                          // Accessing as a bool is likely to undefined
                          // behavior (which triggers a warning from
                          // the sanitizer).  As a char, it may not
                          // be printable.
                          //<< "bool=" << a.second.b
                          //<< " char=" << a.second.c
                          << " int=" << a.second.i
                          << " float=" << a.second.f << " double=" << a.second.d
                          << " void*=" << a.second.p);
    }
  }
}

std::uint32_t TFCSSimulationState::getAuxIndex(const std::string &s) {
  return TFCSSimulationState::fnv1a_32(s.c_str(), s.size());
}

std::uint32_t TFCSSimulationState::getAuxIndex(const char *s) {
  return TFCSSimulationState::fnv1a_32(s, std::strlen(s));
}

void TFCSSimulationState::AddAuxInfoCleanup(
    const TFCSParametrizationBase *para) {
  m_AuxInfoCleanup.insert(para);
}

void TFCSSimulationState::DoAuxInfoCleanup() {
  for (const auto *para : m_AuxInfoCleanup) {
    para->CleanAuxInfo(*this);
  }
}

#pragma GCC diagnostic pop 
