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
  {
  }

  // Constructor with particle container and label
  explicit Event(ParticleContainer container, std::string label)
      : m_ptcl_container(std::move(container))
      , m_label(std::move(label))
  {
  }

  void set_label(std::string label) noexcept { m_label = std::move(label); }

  auto get_label() const -> std::string { return m_label; }

  auto get_container() const -> const ParticleContainer&
  {
    return m_ptcl_container;
  }

private:
  ParticleContainer m_ptcl_container;
  std::string m_label;
};

}  // namespace TestHelpers