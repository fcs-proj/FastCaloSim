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
  using Point = bg::model::point<double, 3, bg::cs::cartesian>;
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

  std::size_t size() const { return m_cells.size(); }
  const Cell& at(std::size_t idx) const { return m_cells.at(idx); }

private:
  bgi::rtree<Value, bgi::quadratic<rtree_quadratic>> m_rtree;
  std::vector<Cell> m_cells;
  CoordinateSystem m_coordinate_system = CoordinateSystem::Undefined;

  static auto determine_coordinate_system(const Cell& cell) -> CoordinateSystem
  {
    if (cell->isXYZ())
      return CoordinateSystem::XYZ;
    if (cell->isEtaPhiR())
      return CoordinateSystem::EtaPhiR;
    if (cell->isEtaPhiZ())
      return CoordinateSystem::EtaPhiZ;
    return CoordinateSystem::Undefined;
  }

  template<typename T>
  auto create_query_point(const T& hit) const -> Point
  {
    switch (m_coordinate_system) {
      case CoordinateSystem::EtaPhiR:
        return Point(hit.eta(), hit.phi(), hit.r());
      case CoordinateSystem::EtaPhiZ:
        return Point(hit.eta(), hit.phi(), hit.z());
      case CoordinateSystem::XYZ:
        return Point(hit.x(), hit.y(), hit.z());
      default:
        throw std::logic_error("Undefined coordinate system");
    }
  }

  static auto normalize_angle(double angle) -> double
  {
    angle = std::fmod(angle + M_PI, 2.0 * M_PI);
    if (angle < 0) {
      angle += 2.0 * M_PI;
    }
    return angle - M_PI;
  }

  void insert(const Cell& cell)
  {
    m_cells.push_back(cell);

    const float half = 0.5;
    std::vector<Box> boxes;
    if (cell->isXYZ()) {
      double xmin = cell->x() - cell->dx() * half;
      double xmax = cell->x() + cell->dx() * half;
      double ymin = cell->y() - cell->dy() * half;
      double ymax = cell->y() + cell->dy() * half;
      double zmin = cell->z() - cell->dz() * half;
      double zmax = cell->z() + cell->dz() * half;
      boxes.emplace_back(Point(xmin, ymin, zmin), Point(xmax, ymax, zmax));
      insert_boxes(cell, boxes);
      return;
    }

    double etamin = cell->eta() - cell->deta() * half;
    double etamax = cell->eta() + cell->deta() * half;
    double phimin = normalize_angle(cell->phi() - cell->dphi() * half);
    double phimax = normalize_angle(cell->phi() + cell->dphi() * half);

    if (cell->isEtaPhiR()) {
      double rmin = cell->r() - cell->dr() * half;
      double rmax = cell->r() + cell->dr() * half;
      if (phimin > phimax) {
        boxes.emplace_back(Point(etamin, -M_PI, rmin),
                           Point(etamax, phimax, rmax));
        boxes.emplace_back(Point(etamin, phimin, rmin),
                           Point(etamax, M_PI, rmax));
      } else {
        boxes.emplace_back(Point(etamin, phimin, rmin),
                           Point(etamax, phimax, rmax));
      }
    } else if (cell->isEtaPhiZ()) {
      double zmin = cell->z() - cell->dz() * half;
      double zmax = cell->z() + cell->dz() * half;
      if (phimin > phimax) {
        boxes.emplace_back(Point(etamin, -M_PI, zmin),
                           Point(etamax, phimax, zmax));
        boxes.emplace_back(Point(etamin, phimin, zmin),
                           Point(etamax, M_PI, zmax));
      } else {
        boxes.emplace_back(Point(etamin, phimin, zmin),
                           Point(etamax, phimax, zmax));
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
