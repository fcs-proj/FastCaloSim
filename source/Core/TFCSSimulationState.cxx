// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#include <cstring>
#include <iostream>
#include <type_traits>
#include <variant>

#include "FastCaloSim/Core/TFCSSimulationState.h"

#include "CLHEP/Random/RandomEngine.h"
#include "FastCaloSim/Definitions/FloatCompare.h"

//=============================================
//======= TFCSSimulationState =========
//=============================================

TFCSSimulationState::TFCSSimulationState(CLHEP::HepRandomEngine* randomEngine)
    : m_randomEngine(randomEngine)
{
  clear();
}

void TFCSSimulationState::clear()
{
  set_SF(1);
  m_Ebin.reset();
  m_Etot = 0;
  m_E.clear();
  m_Efrac.clear();
  m_last_cell.reset();
}

void TFCSSimulationState::deposit(const unsigned long long cell_id, double E)
{
  if (m_last_cell.energy == nullptr || cell_id != m_last_cell.id) {
    m_last_cell.id = cell_id;
    m_last_cell.energy = &m_cells[cell_id];
  }
  *m_last_cell.energy += E;
}

void TFCSSimulationState::Print(Option_t* /*option*/) const
{
  FCS_MSG_INFO("Ebin=" << (m_Ebin ? std::to_string(*m_Ebin) : "none")
                       << " E=" << E() << " #cells=" << m_cells.size());
  for (std::size_t i = 0; i < m_E.size(); ++i) {
    if (!FastCaloSim::exactly_equal(E(i), 0.0)) {
      FCS_MSG_INFO("  E" << i << E(i) << " E" << i << "/E=" << Efrac(i));
    }
  }
  if (!m_AuxInfo.empty()) {
    FCS_MSG_INFO("  AuxInfo has " << m_AuxInfo.size() << " elements");
    for (const auto& a : m_AuxInfo) {
      // Print the value with the type it was stored with
      std::visit(
          [&a, this](const auto& value) -> void
          {
            // Print chars as numbers, they may not be printable
            if constexpr (std::is_same_v<std::decay_t<decltype(value)>, char>) {
              FCS_MSG_INFO("    " << a.first << " : "
                                  << static_cast<int>(value));
            } else {
              FCS_MSG_INFO("    " << a.first << " : " << value);
            }
          },
          a.second);
    }
  }
  if (!m_AuxObjects.empty()) {
    FCS_MSG_INFO("  AuxObjects has " << m_AuxObjects.size() << " elements");
    for (const auto& a : m_AuxObjects) {
      FCS_MSG_INFO("    " << a.first);
    }
  }
}

auto TFCSSimulationState::getAuxIndex(const std::string& s) -> std::uint32_t
{
  return TFCSSimulationState::fnv1a_32(s.c_str(), s.size());
}

auto TFCSSimulationState::getAuxIndex(const char* s) -> std::uint32_t
{
  return TFCSSimulationState::fnv1a_32(s, std::strlen(s));
}
