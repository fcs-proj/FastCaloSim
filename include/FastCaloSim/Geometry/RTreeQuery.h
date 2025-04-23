// Copyright (c) 2025 CERN for the benefit of the FastCaloSim project
#pragma once

#include <memory>
#include <string>

#include <spatialindex/SpatialIndex.h>

#include "FastCaloSim/Geometry/Cell.h"
#include "FastCaloSim/Geometry/RTreeHelpers.h"

/**
 * @brief Class responsible for loading and querying R-trees
 *
 * This class handles loading R-trees from disk with configurable
 * caching for efficient lookup of the nearest cell for a given position.
 */
class RTreeQuery
{
public:
  /**
   * @brief Construct an RTreeQuery with a specific coordinate system
   * @param coordSys The coordinate system used by this layer's cells
   */
  explicit RTreeQuery(RTreeHelpers::CoordinateSystem coordSys);

  ~RTreeQuery() = default;

  /**
   * @brief Load an R-tree from disk with caching
   * @param base_path Path to the R-tree index file
   * @param cache_size Size of the cache in bytes (default: 256KB)
   */
  void load(const std::string& base_path, size_t cache_size = 262144);

  /**
   * @brief Query the nearest cell for a given position
   * @param pos The position to query
   * @return Pointer to the nearest cell
   */
  auto query_point(const Position& pos) const -> uint64_t;

private:
  RTreeHelpers::CoordinateSystem m_coordinate_system;
  std::unique_ptr<SpatialIndex::ISpatialIndex> m_tree;
  std::unique_ptr<SpatialIndex::IStorageManager> m_diskfile;
};
