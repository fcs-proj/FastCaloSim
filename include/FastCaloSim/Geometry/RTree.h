// Copyright (c) 2025 CERN for the benefit of the FastCaloSim project

#pragma once

#include <stdexcept>
#include <vector>

#include <boost/geometry.hpp>
#include <boost/geometry/index/rtree.hpp>

#include "FastCaloSim/Geometry/Cell.h"
namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;
class RTree
{
  using Point = bg::model::point<double, 2, bg::cs::cartesian>;
  using Box = bg::model::box<Point>;
  using Value = std::pair<Box, Cell>;
  // Note: We use rstar as the parameter type, but the tree will be created
  // using the packing algorithm when we pass the full range of values to the
  // constructor See also:
  // https://www.boost.org/doc/libs/1_88_0/boost/geometry/index/rtree.hpp
  using RtreeType = bgi::rtree<Value, bgi::rstar<16>>;

  enum class CoordinateSystem
  {
    XYZ,
    EtaPhiR,
    EtaPhiZ,
    RPhiZ,
    Undefined
  };

public:
  void insert_cell(const Cell& cell)
  {
    CoordinateSystem coordinate_system = determine_coordinate_system(cell);
    if (m_coordinate_system == CoordinateSystem::Undefined) {
      m_coordinate_system = coordinate_system;
    } else if (m_coordinate_system != coordinate_system) {
      throw std::invalid_argument(
          "Tree can't be built with cells with differing coordinate systems!");
    }

    // Create and store boxes for the cell
    prepare_boxes_for_cell(cell);
  }

  // Build the rtree with all stored cells using the packing algorithm
  void build()
  {
    if (m_values.empty()) {
      return;
    }
    // Create a new rtree passing all values at once, which will use the packing
    // algorithm
    m_rtree = std::make_unique<RtreeType>(m_values.begin(), m_values.end());
  }

  ///
  /// @brief Returns the nearest cell (bounding box) for a given hit
  /// Note: The hit must have the same coordinate system as the cells in the
  /// tree
  ///
  auto query_point(const Position& pos) const -> Cell
  {
    if (!m_rtree) {
      throw std::logic_error("Tree not built yet. Call build() first.");
    }
    Point query_point = create_query_point(pos);
    std::vector<Value> result;
    m_rtree->query(bgi::nearest(query_point, 1), std::back_inserter(result));
    if (!result.empty()) {
      return result[0].second;
    }
    throw std::invalid_argument("No cell found for query point");
  }

private:
  std::unique_ptr<RtreeType> m_rtree;
  // Store all cells and values for bulk insertion into the rtree
  std::vector<Value> m_values;
  CoordinateSystem m_coordinate_system = CoordinateSystem::Undefined;

  static auto determine_coordinate_system(const Cell& cell) -> CoordinateSystem
  {
    return cell.isXYZ()    ? CoordinateSystem::XYZ
        : cell.isEtaPhiR() ? CoordinateSystem::EtaPhiR
        : cell.isEtaPhiZ() ? CoordinateSystem::EtaPhiZ
        : cell.isRPhiZ()   ? CoordinateSystem::RPhiZ
                           : CoordinateSystem::Undefined;
  }

  auto create_query_point(const Position& pos) const -> Point
  {
    switch (m_coordinate_system) {
      case CoordinateSystem::EtaPhiR:
      case CoordinateSystem::EtaPhiZ:
      case CoordinateSystem::RPhiZ:
        return {pos.eta(), pos.phi()};
      case CoordinateSystem::XYZ:
        return {pos.x(), pos.y()};
      default:
        throw std::logic_error("Undefined coordinate system");
    }
  }

  void prepare_boxes_for_cell(const Cell& cell)
  {
    const float half = 0.5;
    std::vector<Box> boxes;
    if (cell.isXYZ()) {
      double xmin = cell.x() - cell.dx() * half;
      double xmax = cell.x() + cell.dx() * half;
      double ymin = cell.y() - cell.dy() * half;
      double ymax = cell.y() + cell.dy() * half;
      boxes.emplace_back(Point(xmin, ymin), Point(xmax, ymax));
    } else {
      double etamin = cell.eta() - cell.deta() * half;
      double etamax = cell.eta() + cell.deta() * half;
      double phimin = Cell::norm_angle(cell.phi() - cell.dphi() * half);
      double phimax = Cell::norm_angle(cell.phi() + cell.dphi() * half);

      if (cell.isEtaPhiR() || cell.isEtaPhiZ() || cell.isRPhiZ()) {
        if (phimin > phimax) {
          boxes.emplace_back(Point(etamin, -M_PI), Point(etamax, phimax));
          boxes.emplace_back(Point(etamin, phimin), Point(etamax, M_PI));
        } else {
          boxes.emplace_back(Point(etamin, phimin), Point(etamax, phimax));
        }
      }
    }

    // Store all boxes with their associated cell for later bulk insertion
    for (const auto& box : boxes) {
      m_values.emplace_back(box, cell);
    }
  }
};
