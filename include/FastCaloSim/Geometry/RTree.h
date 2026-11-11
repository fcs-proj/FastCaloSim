// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project
#pragma once

#include <cstdint>
#include <memory>
#include <vector>

#include <FastCaloSim/FastCaloSim_export.h>

#include "FastCaloSim/Geometry/Cell.h"
#include "FastCaloSim/Geometry/RTreeHelpers.h"

/**
 * @brief Finds the cell of a layer that a position belongs to
 *
 * In the lookup coordinates of its layer, (x, y) or (eta, phi), every cell
 * is a rectangle. The rectangles of all cells of the layer are held in an
 * R-tree (Boost.Geometry), which is built once, kept in memory and shared
 * read-only by all threads.
 */
class FASTCALOSIM_EXPORT RTree
{
public:
  /**
   * @brief Build the R-tree of a layer
   * @param coord_sys The coordinate system of the layer's cells
   * @param cells The cells of the layer
   */
  RTree(RTreeHelpers::CoordinateSystem coord_sys,
        const std::vector<Cell>& cells);

  ~RTree();

  /**
   * @brief Find the cell for a given position
   * @param pos The position to query
   * @return ID of the cell whose rectangle contains the position.
   *
   * The rectangles of neighboring cells can overlap slightly. Where several
   * contain the position, the cell with the nearest rectangle center is
   * returned (the smaller ID on an exact tie). Where none contains it (in a
   * gap between cells or outside the layer), the cell with the nearest
   * rectangle is returned.
   */
  auto query_point(const Position& pos) const -> uint64_t;

private:
  // Coordinate system for this layer
  RTreeHelpers::CoordinateSystem m_coordinate_system;

  // The Boost.Geometry R-tree. It is defined in the source file, so that the
  // installed headers do not depend on Boost.
  struct Impl;
  std::unique_ptr<const Impl> m_impl;
};
