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
  static constexpr int rtree_quadratic = 16;

  enum class CoordinateSystem
  {
    XYZ,
    EtaPhiR,
    EtaPhiZ,
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

    insert(cell);
  }

  template<typename T>
  auto query_point(const T& hit) const -> Cell
  {
    Point query_point = create_query_point(hit);

    std::vector<Value> result;
    m_rtree.query(bgi::nearest(query_point, 1), std::back_inserter(result));

    if (!result.empty()) {
      return result[0].second;
    }

    throw std::invalid_argument("No cell found for query point");
  }

  auto size() const -> std::size_t { return m_cells.size(); }
  auto at(std::size_t idx) const -> const Cell& { return m_cells.at(idx); }

private:
  bgi::rtree<Value, bgi::quadratic<rtree_quadratic>> m_rtree;
  std::vector<Cell> m_cells;
  CoordinateSystem m_coordinate_system = CoordinateSystem::Undefined;

  static auto determine_coordinate_system(const Cell& cell) -> CoordinateSystem
  {
    return cell.isXYZ()    ? CoordinateSystem::XYZ
        : cell.isEtaPhiR() ? CoordinateSystem::EtaPhiR
        : cell.isEtaPhiZ() ? CoordinateSystem::EtaPhiZ
                           : CoordinateSystem::Undefined;
  }

  template<typename T>
  auto create_query_point(const T& hit) const -> Point
  {
    switch (m_coordinate_system) {
      case CoordinateSystem::EtaPhiR:
      case CoordinateSystem::EtaPhiZ:
        return {hit.eta(), hit.phi()};
      case CoordinateSystem::XYZ:
        return {hit.x(), hit.y()};
      default:
        throw std::logic_error("Undefined coordinate system");
    }
  }

  void insert(const Cell& cell)
  {
    m_cells.push_back(cell);

    const float half = 0.5;
    std::vector<Box> boxes;
    if (cell.isXYZ()) {
      double xmin = cell.x() - cell.dx() * half;
      double xmax = cell.x() + cell.dx() * half;
      double ymin = cell.y() - cell.dy() * half;
      double ymax = cell.y() + cell.dy() * half;

      boxes.emplace_back(Point(xmin, ymin), Point(xmax, ymax));
      insert_boxes(cell, boxes);
      return;
    }

    double etamin = cell.eta() - cell.deta() * half;
    double etamax = cell.eta() + cell.deta() * half;
    double phimin = Cell::norm_angle(cell.phi() - cell.dphi() * half);
    double phimax = Cell::norm_angle(cell.phi() + cell.dphi() * half);

    if (cell.isEtaPhiR()) {
      if (phimin > phimax) {
        boxes.emplace_back(Point(etamin, -M_PI), Point(etamax, phimax));
        boxes.emplace_back(Point(etamin, phimin), Point(etamax, M_PI));
      } else {
        boxes.emplace_back(Point(etamin, phimin), Point(etamax, phimax));
      }
    } else if (cell.isEtaPhiZ()) {
      if (phimin > phimax) {
        boxes.emplace_back(Point(etamin, -M_PI), Point(etamax, phimax));
        boxes.emplace_back(Point(etamin, phimin), Point(etamax, M_PI));
      } else {
        boxes.emplace_back(Point(etamin, phimin), Point(etamax, phimax));
      }
    }
    insert_boxes(cell, boxes);
  }

  void insert_boxes(const Cell& cell, const std::vector<Box>& boxes)
  {
    for (const auto& box : boxes) {
      m_rtree.insert(std::make_pair(box, cell));
    }
  }
};
