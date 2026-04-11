#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-equal"

// Copyright (c) 2024 CERN for the benefit of the FastCaloSim project

#include <cstring>
#include <iostream>

#include "FastCaloSim/Core/TFCSSimulationState.h"

#include "CLHEP/Random/RandomEngine.h"
#include "FastCaloSim/Core/TFCSParametrizationBase.h"

//=============================================
//======= TFCSSimulationState =========
//=============================================
using namespace FastCaloSim::Core;

TFCSSimulationState::TFCSSimulationState(CLHEP::HepRandomEngine* randomEngine)
    : m_randomEngine(randomEngine)
{
  clear();
}

void TFCSSimulationState::clear()
{
  set_SF(1);
  m_Ebin = -1;
  m_Etot = 0;
  m_E.clear();
  m_Efrac.clear();
}

void TFCSSimulationState::deposit(const unsigned long long cell_id, float E)
{
  m_cells[cell_id] += E;
}

void TFCSSimulationState::Print(Option_t*) const
{
  MSG_INFO("Ebin=" << m_Ebin << " E=" << E() << " #cells=" << m_cells.size());
  for (int i = 0; i < m_E.size(); ++i)
    if (E(i) != 0) {
      MSG_INFO("  E" << i << E(i) << " E" << i << "/E=" << Efrac(i));
    }
  if (!m_AuxInfo.empty()) {
    MSG_INFO("  AuxInfo has " << m_AuxInfo.size() << " elements");
    for (const auto& a : m_AuxInfo) {
      MSG_INFO("    " << a.first
                      << " : "
                      // Don't print as char/bool.
                      // Accessing as a bool is likely to undefined
                      // behavior (which triggers a warning from
                      // the sanitizer).  As a char, it may not
                      // be printable.
                      //<< "bool=" << a.second.b
                      //<< " char=" << a.second.c
                      << " int=" << a.second.i << " float=" << a.second.f
                      << " double=" << a.second.d << " void*=" << a.second.p);
    }
  }
}

std::uint32_t TFCSSimulationState::getAuxIndex(const std::string& s)
{
  return TFCSSimulationState::fnv1a_32(s.c_str(), s.size());
}

std::uint32_t TFCSSimulationState::getAuxIndex(const char* s)
{
  return TFCSSimulationState::fnv1a_32(s, std::strlen(s));
}

void TFCSSimulationState::AddAuxInfoCleanup(const TFCSParametrizationBase* para)
{
  m_AuxInfoCleanup.insert(para);
}

void TFCSSimulationState::DoAuxInfoCleanup()
{
  for (const auto* para : m_AuxInfoCleanup) {
    para->CleanAuxInfo(*this);
  }
}

#pragma GCC diagnostic pop
