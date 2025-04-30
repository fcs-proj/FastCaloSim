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
  try {
    std::string filename = base_path;
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

// Internal reusable visitor (anonymous namespace)
namespace
{
class CellVisitor : public SpatialIndex::IVisitor
{
public:
  CellVisitor(bool match_id = false, unsigned long long id = 0)
      : m_match_id(match_id)
      , m_target_id(id)
  {
  }

  Cell result;
  bool found = false;

  void visitData(const SpatialIndex::IData& d) override
  {
    if (m_match_id && d.getIdentifier() != m_target_id)
      return;

    uint8_t* buffer = nullptr;
    uint32_t len = 0;
    d.getData(len, &buffer);

    if (len == sizeof(Cell)) {
      std::memcpy(&result, buffer, sizeof(Cell));
      found = true;
    }
  }

  void visitNode(const SpatialIndex::INode&) override {}
  void visitData(std::vector<const SpatialIndex::IData*>&) override {}

private:
  bool m_match_id;
  unsigned long long m_target_id;
};
}  // namespace

/**
 * @brief Retrieve the cell closest to a given position in 2D (either XY or
 * eta-phi).
 *
 * This function performs a nearest neighbor query on the R-tree using the
 * provided position. It returns the first matching cell whose centroid is
 * closest to the input coordinates according to the coordinate system
 * configured (XY or eta-phi).
 *
 * The search is done using `libspatialindex`'s `nearestNeighborQuery()` method,
 * which is spatially indexed and highly efficient — suitable for fast, repeated
 * use in simulation
 *
 * @note If the coordinate system is `CoordinateSystem::XYZ`, the query is
 * performed in the XY plane. If it is `CoordinateSystem::ETAPHI`, then the
 * query is performed in eta-phi space.
 *
 * @warning This function assumes that the R-tree has been correctly populated
 * with 2D bounding boxes using matching coordinates (e.g. XY or eta-phi).
 * Querying with mismatched coordinate conventions may return incorrect results.
 *
 * @param pos The position to query (must provide x/y or eta/phi depending on
 * the coordinate mode).
 * @return The closest Cell to the specified point, or a default-constructed
 * Cell if none found.
 */
auto RTreeQuery::query_point(const Position& pos) const -> Cell
{
  if (!m_tree) {
    throw std::logic_error("Tree not loaded yet. Call load() first.");
  }

  double coords[2];
  if (m_coordinate_system == RTreeHelpers::CoordinateSystem::XYZ) {
    coords[0] = pos.x();
    coords[1] = pos.y();
  } else {
    coords[0] = pos.eta();
    coords[1] = pos.phi();
  }

  SpatialIndex::Point query(coords, 2);
  CellVisitor visitor;

  try {
    m_tree->nearestNeighborQuery(1, query, visitor);
    if (visitor.found) {
      return visitor.result;
    }
  } catch (Tools::Exception& e) {
    throw std::runtime_error(std::string("Error querying RTree: ") + e.what());
  }

  return Cell();  // Invalid/default Cell
}

/**
 * @brief Retrieve a cell by its unique identifier (ID).
 *
 * This method performs a full scan of the R-tree using a bounding region
 * that encompasses the entire spatial index and filters entries by ID.
 *
 * @warning This function is **slow**: it requires a full traversal of the
 * R-tree, because libspatialindex does not provide a direct index or hash map
 * for fast lookup by identifier.
 *
 * @note This method is intended for debugging, visualization, or rare use cases
 * where position-based queries are not feasible. It should **not** be used in
 * performance-sensitive contexts such as full simulation, where thousands to
 * millions of cells may need to be queried. In such cases, positional queries
 * (e.g., `query_point()`) or precomputed data structures (e.g.,
 * `unordered_map<id, Cell>`) should be used instead.
 *
 * @param id The unique cell identifier (same as what was stored as the R-tree
 * entry ID).
 * @return The matching Cell, or a default/invalid Cell if not found.
 */
auto RTreeQuery::get_cell_by_id(unsigned long long id) const -> Cell
{
  if (!m_tree) {
    throw std::logic_error("Tree not loaded yet. Call load() first.");
  }

  CellVisitor visitor(true, id);

  double low[2] = {-1e10, -1e10};
  double high[2] = {1e10, 1e10};
  SpatialIndex::Region universe(low, high, 2);

  try {
    m_tree->intersectsWithQuery(universe, visitor);
    if (visitor.found) {
      return visitor.result;
    }
  } catch (Tools::Exception& e) {
    throw std::runtime_error(std::string("Error during ID lookup: ")
                             + e.what());
  }

  return Cell();  // Not found
}
