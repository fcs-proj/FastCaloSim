// Copyright (c) 2024 CERN for the benefit of the FastCaloSim project

#pragma once

#include <utility>

#include "ParticleContainer.h"

namespace TestHelpers
{

class Event
{
public:
  // Default constructor
  Event() = default;
  // Constructor with particle container
  explicit Event(ParticleContainer container)
      : m_ptcl_container(std::move(container))
      , m_event_id(-1)
      , m_print_string("Event")
  {
  }

  // Constructor with particle container and label
  explicit Event(ParticleContainer container, std::string label)
      : m_ptcl_container(std::move(container))
      , m_label(std::move(label))
      , m_event_id(-1)
  {
  }

  // Event label setter and getter
  inline auto set_label(std::string label) -> void
  {
    m_label = std::move(label);
  }
  inline auto get_label() const -> std::string { return m_label; }

  // Print string setter and getter
  inline auto set_print_string(std::string print_string) -> void
  {
    m_print_string = std::move(print_string);
  }
  inline auto get_print_string() const -> std::string { return m_print_string; }

  // Event ID setter and getter
  inline auto set_event_id(int event_id) -> void { m_event_id = event_id; }
  inline auto get_event_id() const -> int { return m_event_id; }

  // Retrieve particle container associated to event
  inline auto get_container() const -> const ParticleContainer&
  {
    return m_ptcl_container;
  }

  // Set the param file name used to simulate the event
  inline auto set_param_file_key(const std::string& param_file_key) -> void
  {
    m_param_file_key = param_file_key;
  }

  // Get the param file name used to simulate the event
  inline auto get_param_file_key() const -> std::string
  {
    return m_param_file_key;
  }

private:
  // @brief Container holding the particles
  ParticleContainer m_ptcl_container;
  // @brief Label of the event
  std::string m_label;
  // @brief Identifier of the event
  int m_event_id;
  // @brief String used by gtest to name event
  std::string m_print_string;
  // @brief Key of the param file used to simulate the event
  std::string m_param_file_key;

  // Friend method to let google test know how to label a test parametrized by
  // If event_id is set, the print string will be the label + "_" + event_id
  friend void PrintTo(const Event& evt, std::ostream* os)
  {
    if (evt.m_event_id == -1) {
      *os << evt.m_print_string;
    } else {
      *os << evt.m_print_string + "_" + std::to_string(evt.m_event_id);
    }
  }
};

}  // namespace TestHelpers
