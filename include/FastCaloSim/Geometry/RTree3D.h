#include <vector>

#include <boost/geometry.hpp>
#include <boost/geometry/index/rtree.hpp>

#include "FastCaloSim/Geometry/Cell.h"

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

using Point3D = bg::model::point<double, 3, bg::cs::cartesian>;
using Box = bg::model::box<Point3D>;
using Value = std::pair<Box, Cell>;

class RTree3D
{
  enum class CoordinateSystem
  {
    Cartesian,
    Cylindrical,
    ECCylindrical,
    Undefined
  };

private:
  bgi::rtree<Value, bgi::quadratic<16>> m_rtree;
  std::vector<Cell> m_cells;
  CoordinateSystem m_coordinate_system = CoordinateSystem::Undefined;

  static auto normalize_angle(double angle) -> double
  {
    angle = std::fmod(angle + M_PI, 2.0 * M_PI);
    if (angle < 0) {
      angle += 2.0 * M_PI;
    }
    return angle - M_PI;
  }

  // Function to split bounding boxes that wrap around -π to π
  // the logic here needs to be revisited but I think we are on the right track
  void insert_wrapped_bounding_boxes(const Cell& cell,
                                     double eta_min,
                                     double eta_max,
                                     double phi_min,
                                     double phi_max,
                                     double r_min,
                                     double r_max)
  {
    // Normalize phi_min and phi_max
    phi_min = normalize_angle(phi_min);
    phi_max = normalize_angle(phi_max);

    if (phi_min > phi_max)
    {  // needs to be revisited, but I think this looks OK
      // Bounding box wraps around, split into two boxes
      Box box1 = {Point3D {eta_min, -M_PI, r_min},
                  Point3D {eta_max, phi_max, r_max}};
      Box box2 = {Point3D {eta_min, phi_min, r_min},
                  Point3D {eta_max, M_PI, r_max}};
      m_rtree.insert(std::make_pair(box1, cell));
      m_rtree.insert(std::make_pair(box2, cell));
    } else {
      // Single bounding box
      Box box = {Point3D {eta_min, phi_min, r_min},
                 Point3D {eta_max, phi_max, r_max}};
      m_rtree.insert(std::make_pair(box, cell));
    }
  }

  auto get_bounding_box(const Cell& cell) -> void
  {
    if (cell.isCylindrical || cell.isECCylindrical) {
      double eta_min = cell.eta - cell.deta / 2.0;
      double eta_max = cell.eta + cell.deta / 2.0;
      double phi_min = cell.phi - cell.dphi / 2.0;
      double phi_max = cell.phi + cell.dphi / 2.0;
      double r_min =
          cell.isCylindrical ? cell.r - cell.dr / 2.0 : cell.z - cell.dz / 2.0;
      double r_max =
          cell.isCylindrical ? cell.r + cell.dr / 2.0 : cell.z + cell.dz / 2.0;

      insert_wrapped_bounding_boxes(
          cell, eta_min, eta_max, phi_min, phi_max, r_min, r_max);
    } else if (cell.isCartesian) {
      double x_min = cell.x - cell.dx / 2.0;
      double x_max = cell.x + cell.dx / 2.0;
      double y_min = cell.y - cell.dy / 2.0;
      double y_max = cell.y + cell.dy / 2.0;
      double z_min = cell.z - cell.dz / 2.0;
      double z_max = cell.z + cell.dz / 2.0;
      Box box = {Point3D {x_min, y_min, z_min}, Point3D {x_max, y_max, z_max}};
      m_rtree.insert(std::make_pair(box, cell));
    } else {
      throw std::invalid_argument("Unknown cell coordinate system");
    }
  }

public:
  void insert_cell(const Cell& cell)
  {
    CoordinateSystem coordinate_system = cell.isCartesian
        ? CoordinateSystem::Cartesian
        : cell.isCylindrical ? CoordinateSystem::Cylindrical
                             : CoordinateSystem::ECCylindrical;

    if (m_coordinate_system == CoordinateSystem::Undefined) {
      m_coordinate_system = coordinate_system;
    } else if (m_coordinate_system != coordinate_system) {
      std::cout << " Cell trying to be added has coordinate system "
                << static_cast<int>(coordinate_system)
                << " but tree has coordinate system "
                << static_cast<int>(m_coordinate_system) << '\n';

      throw std::invalid_argument(
          "Tree can't be built with cells with differing coordinate systems!");
    }

    get_bounding_box(cell);
    m_cells.push_back(cell);
  }

  // Converts the query into the coordinate system the tree was built with
  auto convert(double x, double y, double z) -> Point3D
  {
    Point3D query_point;
    if (m_coordinate_system == CoordinateSystem::Cartesian) {
      query_point = Point3D {x, y, z};
    } else if (m_coordinate_system == CoordinateSystem::Cylindrical) {
      // Convert to eta, phi, r
      double r = std::sqrt(x * x + y * y);
      double phi = normalize_angle(std::atan2(y, x));
      double theta = std::atan2(r, z);
      double eta = -std::log(std::tan(theta / 2.0));
      query_point = Point3D {eta, phi, r};
    } else if (m_coordinate_system == CoordinateSystem::ECCylindrical) {
      // Convert to eta, phi, z
      double r = std::sqrt(x * x + y * y);
      double phi = normalize_angle(std::atan2(y, x));
      double theta = std::atan2(r, z);
      double eta = -std::log(std::tan(theta / 2.0));
      query_point = Point3D {eta, phi, z};
    } else {
      throw std::invalid_argument("Unknown cell coordinate system");
    }
    return query_point;
  }

  auto query_point(double x, double y, double z) -> Cell
  {
    Point3D query_point = convert(x, y, z);

    std::vector<Value> result;
    m_rtree.query(bgi::nearest(query_point, 1), std::back_inserter(result));

    if (!result.empty()) {
      return result[0].second;
    }

    // Return an invalid cell to indicate no cell was found
    return {-1,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            false,
            false,
            false,
            false,
            false};
  }

  auto size() const -> std::size_t { return m_cells.size(); }
  auto at(std::size_t idx) const -> const Cell& { return m_cells.at(idx); }
};
