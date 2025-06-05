// Copyright (c) 2025 CERN for the benefit of the FastCaloSim project
#pragma once

#include <array>
#include <memory>
#include <string>
#include <vector>

#include <spatialindex/SpatialIndex.h>

#include "FastCaloSim/Geometry/Cell.h"
#include "FastCaloSim/Geometry/RTreeHelpers.h"
namespace FastCaloSim::Geometry
{

/**
 * @brief Class responsible for building and persisting R-trees
 *
 * This class handles collecting cells and bulk loading them into an R-tree
 * for a specific detector layer with a known coordinate system.
 */
class RTreeBuilder
{
public:
  /**
   * @brief Construct an RTreeBuilder with a specific coordinate system
   * @param coordSys The coordinate system for this layer's cells
   */
  explicit RTreeBuilder(RTreeHelpers::CoordinateSystem coordSys);

  ~RTreeBuilder() = default;

  /**
   * @brief Add a cell for inclusion in the R-tree
   * @param cell Pointer to the cell to add
   */
  void add_cell(const Cell* cell);

  /**
   * @brief Build the R-tree using bulk loading and persist to disk
   * @param output_path Path where the R-tree index will be stored
   */
  void build(const std::string& output_path);

private:
  RTreeHelpers::CoordinateSystem m_coordinate_system;
  // Store bounding boxes and their associated cell IDs
  std::vector<std::pair<std::array<double, 4>, uint64_t>> m_boxes;
};
}  // namespace FastCaloSim::Geometry
