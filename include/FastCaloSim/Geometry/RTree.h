// Copyright (c) 2025 CERN for the benefit of the FastCaloSim project
#pragma once

#include <memory>
#include <string>

#include "FastCaloSim/Geometry/Cell.h"
#include "FastCaloSim/Geometry/RTreeBuilder.h"
#include "FastCaloSim/Geometry/RTreeHelpers.h"
#include "FastCaloSim/Geometry/RTreeQuery.h"

/**
 * @brief Main R-tree class providing a unified interface
 *
 * This class delegates to specialized components for building,
 * persisting, loading, and querying R-trees. It provides a simple
 * interface for working with R-trees at the layer level with a
 * known coordinate system.
 *
 * The class supports:
 * 1. Bulk loading of R-tree for better performance
 * 2. Memory persistence capability
 * 3. Efficient caching for fast lookups
 */
class RTree
{
public:
  /**
   * @brief Construct an RTree with a specific coordinate system
   * @param coordSys The coordinate system for this layer's cells
   */
  explicit RTree(RTreeHelpers::CoordinateSystem coordSys);

  ~RTree() = default;

  /**
   * @brief Insert a cell into the collection for bulk loading
   * @param cell Pointer to the cell to insert
   */
  void insert_cell(const Cell& cell);

  /**
   * @brief Build the RTree using bulk loading and store to disk
   * @param base_path Path to store the R-tree
   */
  void build(const std::string& base_path);

  /**
   * @brief Load a previously saved RTree from file with caching
   * @param base_path Path to the input file
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
  // Coordinate system for this layer
  RTreeHelpers::CoordinateSystem m_coordinate_system;

  // Separate components for building and querying
  std::unique_ptr<RTreeBuilder> m_builder;
  std::unique_ptr<RTreeQuery> m_query;
};
