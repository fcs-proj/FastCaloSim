// Copyright (c) 2025 CERN for the benefit of the FastCaloSim project
#pragma once

#include <cstddef>
#include <memory>
#include <string>

#include <spatialindex/SpatialIndex.h>
#include <tbb/enumerable_thread_specific.h>

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
  std::string m_base_path;
  std::size_t m_cache_size {262144};

  // Per-thread view of the shared, read-only on-disk index files.
  // Declaration order matters: tree depends on buffer depends on diskfile,
  // so they must destruct in reverse (tree, buffer, diskfile) — which this
  // member order gives.
  struct TreeHandle
  {
    std::unique_ptr<SpatialIndex::IStorageManager> diskfile;
    std::unique_ptr<SpatialIndex::IStorageManager> buffer;
    std::unique_ptr<SpatialIndex::ISpatialIndex> tree;
  };
  mutable tbb::enumerable_thread_specific<TreeHandle> m_perThread;

  // Lazily initialize and return this thread's handle.
  TreeHandle& localTree() const;
};
