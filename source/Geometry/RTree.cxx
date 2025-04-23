// Copyright (c) 2025 CERN for the benefit of the FastCaloSim project
#include <stdexcept>

#include "FastCaloSim/Geometry/RTree.h"

RTree::RTree(RTreeHelpers::CoordinateSystem coordSys)
    : m_coordinate_system(coordSys)
    , m_builder(std::make_unique<RTreeBuilder>(coordSys))
{
}

void RTree::insert_cell(const Cell* cell)
{
  if (!m_builder) {
    m_builder = std::make_unique<RTreeBuilder>(m_coordinate_system);
  }

  m_builder->add_cell(cell);
}

void RTree::build(const std::string& base_path)
{
  if (!m_builder) {
    throw std::logic_error(
        "No cells added to the RTree. Call insert_cell() first.");
  }

  m_builder->build(base_path);

  // Release builder resources after building
  m_builder.reset();
}

void RTree::load(const std::string& base_path, size_t cache_size)
{
  if (!m_query) {
    m_query = std::make_unique<RTreeQuery>(m_coordinate_system);
  }

  m_query->load(base_path, cache_size);
}

auto RTree::query_point(const Position& pos) const -> uint64_t
{
  if (!m_query) {
    throw std::logic_error("Tree not loaded yet. Call load() first.");
  }

  return m_query->query_point(pos);
}
