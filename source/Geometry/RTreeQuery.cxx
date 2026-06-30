// Copyright (c) 2025 CERN for the benefit of the FastCaloSim project
#include <cstring>
#include <stdexcept>

#include "FastCaloSim/Geometry/RTreeQuery.h"

RTreeQuery::RTreeQuery(RTreeHelpers::CoordinateSystem coordSys)
    : m_coordinate_system(coordSys)
{
}

void RTreeQuery::load(const std::string& base_path, size_t cache_size)
{
  m_base_path = base_path;
  m_cache_size = cache_size;
  localTree();  // validate now; throws on bad path as before
}

RTreeQuery::TreeHandle& RTreeQuery::localTree() const
{
  TreeHandle& h = m_perThread.local();
  if (!h.tree) {
    try {
      // Copy: loadDiskStorageManager's API takes a non-const std::string&,
      // though it does not conceptually modify the path.
      std::string filename = m_base_path;
      h.diskfile.reset(
          SpatialIndex::StorageManager::loadDiskStorageManager(filename));
      h.buffer.reset(
          SpatialIndex::StorageManager::createNewRandomEvictionsBuffer(
              *h.diskfile, m_cache_size, false));
      SpatialIndex::id_type indexId = 1;
      h.tree.reset(SpatialIndex::RTree::loadRTree(*h.buffer, indexId));
    } catch (Tools::Exception& e) {
      throw std::runtime_error(std::string("Error loading RTree: ") + e.what());
    }
  }
  return h;
}

auto RTreeQuery::query_point(const Position& pos) const -> uint64_t
{
  auto& tree = *localTree().tree;

  // Convert position to query point
  double coords[2];
  if (m_coordinate_system == RTreeHelpers::CoordinateSystem::XYZ) {
    coords[0] = pos.x();
    coords[1] = pos.y();
  } else {
    // For eta-phi based coordinate systems
    coords[0] = pos.eta();
    coords[1] = pos.phi();
  }

  // Create query point
  SpatialIndex::Point query(coords, 2);

  // Setup nearest neighbor visitor
  class Visitor : public SpatialIndex::IVisitor
  {
  public:
    uint64_t result = 0;
    bool found = false;

    void visitData(const SpatialIndex::IData& d) override
    {
      result = d.getIdentifier();
      found = true;
    }

    void visitNode(const SpatialIndex::INode&) override {}
    void visitData(std::vector<const SpatialIndex::IData*>&) override {}
  };

  Visitor visitor;

  try {
    // Execute nearest neighbor query
    tree.nearestNeighborQuery(1, query, visitor);

    // Check if we got a result
    if (visitor.found) {
      return visitor.result;  // Return the cell ID directly
    }
  } catch (Tools::Exception& e) {
    throw std::runtime_error(std::string("Error querying RTree: ") + e.what());
  }

  throw std::invalid_argument("No cell found for query point");
}
