// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project
#include <algorithm>
#include <iterator>
#include <stdexcept>
#include <utility>

#include "FastCaloSim/Geometry/RTree.h"

#include <boost/container/small_vector.hpp>
#include <boost/geometry.hpp>
#include <boost/geometry/index/rtree.hpp>

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

namespace
{
using Point = bg::model::point<double, 2, bg::cs::cartesian>;
using Box = bg::model::box<Point>;
// A cell rectangle and the ID of its cell
using Entry = std::pair<Box, uint64_t>;
using Tree = bgi::rtree<Entry, bgi::rstar<16>>;
}  // namespace

struct RTree::Impl
{
  Tree tree;
};

RTree::RTree(RTreeHelpers::CoordinateSystem coord_sys,
             const std::vector<Cell>& cells)
    : m_coordinate_system(coord_sys)
{
  const bool is_xyz = coord_sys == RTreeHelpers::CoordinateSystem::XYZ;

  std::vector<Entry> entries;
  for (const Cell& cell : cells) {
    // One rectangle per cell, or two if the cell crosses phi = +-pi
    const auto rectangles = is_xyz
        ? std::vector {RTreeHelpers::build_xy_box(
              cell.x(), cell.y(), cell.dx(), cell.dy())}
        : RTreeHelpers::build_eta_phi_boxes(
              cell.eta(), cell.phi(), cell.deta(), cell.dphi());

    for (const auto& r : rectangles) {
      entries.emplace_back(Box(Point(r[0], r[1]), Point(r[2], r[3])),
                           cell.id());
    }
  }

  // Constructing the tree from all rectangles at once packs it (bulk loading)
  m_impl =
      std::make_unique<const Impl>(Impl {Tree(entries.begin(), entries.end())});
}

RTree::~RTree() = default;

auto RTree::query_point(const Position& pos) const -> uint64_t
{
  const Point point = m_coordinate_system == RTreeHelpers::CoordinateSystem::XYZ
      ? Point(pos.x(), pos.y())
      : Point(pos.eta(), pos.phi());

  // The cells whose rectangle contains the point: normally exactly one, but
  // the rectangles of neighboring cells can overlap slightly
  boost::container::small_vector<Entry, 4> candidates;
  m_impl->tree.query(bgi::intersects(point), std::back_inserter(candidates));

  // In a gap between cells or outside the layer: the nearest rectangle
  if (candidates.empty()) {
    m_impl->tree.query(bgi::nearest(point, 1), std::back_inserter(candidates));
  }
  if (candidates.empty()) {
    throw std::invalid_argument("No cell found for query point");
  }

  // Of several candidates, take the cell whose rectangle center is nearest
  // (and the smaller cell ID on an exact tie)
  const auto center_distance_and_id = [&point](const Entry& entry)
  {
    const auto center = bg::return_centroid<Point>(entry.first);
    return std::pair(bg::comparable_distance(point, center), entry.second);
  };
  const auto best = std::min_element(
      candidates.begin(),
      candidates.end(),
      [&](const Entry& a, const Entry& b)
      { return center_distance_and_id(a) < center_distance_and_id(b); });
  return best->second;
}
