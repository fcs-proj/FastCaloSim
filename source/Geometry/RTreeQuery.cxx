// Copyright (c) 2025 CERN for the benefit of the FastCaloSim project
#include <cstring>
#include <stdexcept>

#include "FastCaloSim/Geometry/RTreeQuery.h"

RTreeQuery::RTreeQuery(RTreeHelpers::CoordinateSystem coordSys)
    : m_coordinate_system(coordSys)
{
}

void RTreeQuery::load(const std::string& index_path, size_t cache_size)
{
  try {
    std::string filename = index_path;
    SpatialIndex::IStorageManager* diskfile =
        SpatialIndex::StorageManager::loadDiskStorageManager(filename);

    SpatialIndex::IStorageManager* cache =
        SpatialIndex::StorageManager::createNewRandomEvictionsBuffer(
            *diskfile, cache_size, false);

    SpatialIndex::id_type indexId = 1;
    m_tree.reset(SpatialIndex::RTree::loadRTree(*cache, indexId));
    m_diskfile.reset(cache);
  } catch (Tools::Exception& e) {
    throw std::runtime_error(std::string("Error loading RTree: ") + e.what());
  }
}

auto RTreeQuery::query_point(const Position& pos) const -> uint64_t
{
  if (!m_tree) {
    throw std::logic_error("Tree not loaded yet. Call load() first.");
  }

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
      uint32_t len = 0;
      uint8_t* data = nullptr;
      d.getData(len, &data);
      std::memcpy(&result, data, sizeof(uint64_t));
      found = true;
    }

    void visitNode(const SpatialIndex::INode&) override {}
    void visitData(std::vector<const SpatialIndex::IData*>&) override {}
  };

  Visitor visitor;

  try {
    // Execute nearest neighbor query
    m_tree->nearestNeighborQuery(1, query, visitor);

    // Check if we got a result
    if (visitor.found) {
      return visitor.result;  // Return the cell ID directly
    }
  } catch (Tools::Exception& e) {
    throw std::runtime_error(std::string("Error querying RTree: ") + e.what());
  }

  throw std::invalid_argument("No cell found for query point");
}
