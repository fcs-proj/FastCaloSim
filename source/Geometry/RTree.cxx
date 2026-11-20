// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iterator>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

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
// A cell rectangle and the position of its cell in the constructor's cells
using Entry = std::pair<Box, uint64_t>;
using Tree = bgi::rtree<Entry, bgi::rstar<16>>;

/**
 * @brief Uniform grid over the rectangles of a layer
 *
 * Every rectangle is listed in each bin it overlaps, so the rectangles that
 * contain a point are all listed in the point's bin. The bin of a coordinate
 * is a non-decreasing function of it, which guarantees this also for points
 * on a rectangle edge. Looking up the bin and testing its few rectangles is
 * much faster than descending the R-tree.
 */
class Grid
{
public:
  explicit Grid(const std::vector<Entry>& entries)
  {
    if (entries.empty()) {
      return;
    }

    m_x0 = m_y0 = std::numeric_limits<double>::max();
    double x1 = std::numeric_limits<double>::lowest();
    double y1 = std::numeric_limits<double>::lowest();
    std::vector<double> widths, heights;
    widths.reserve(entries.size());
    heights.reserve(entries.size());
    for (const auto& [box, pos] : entries) {
      m_x0 = std::min(m_x0, box.min_corner().get<0>());
      m_y0 = std::min(m_y0, box.min_corner().get<1>());
      x1 = std::max(x1, box.max_corner().get<0>());
      y1 = std::max(y1, box.max_corner().get<1>());
      widths.push_back(box.max_corner().get<0>() - box.min_corner().get<0>());
      heights.push_back(box.max_corner().get<1>() - box.min_corner().get<1>());
    }
    m_x1 = x1;
    m_y1 = y1;

    // Bins about the size of a typical rectangle, but never more bins than a
    // few per rectangle
    const auto median = [](std::vector<double>& v)
    {
      std::nth_element(v.begin(), v.begin() + v.size() / 2, v.end());
      return v[v.size() / 2];
    };
    const double width = median(widths);
    const double height = median(heights);
    const auto n_bins = [](double extent, double size) -> std::size_t
    {
      constexpr std::size_t max_n = 1U << 15U;
      if (!(extent > 0) || !(size > 0)) {
        return 1;
      }
      // Compare before converting: extent / size can exceed any integer
      const double n = std::ceil(extent / size);
      if (!(n < static_cast<double>(max_n))) {
        return max_n;
      }
      return std::max<std::size_t>(1, static_cast<std::size_t>(n));
    };
    std::size_t nx = n_bins(m_x1 - m_x0, width);
    std::size_t ny = n_bins(m_y1 - m_y0, height);
    const std::size_t max_bins = 4 * entries.size() + 16;
    while (nx * ny > max_bins) {
      nx = std::max<std::size_t>(1, (nx + 1) / 2);
      ny = std::max<std::size_t>(1, (ny + 1) / 2);
    }

    // Rectangles much larger than a bin would be listed in many bins: coarsen
    // the grid until the listing stays small
    while (true) {
      set_bins(nx, ny);
      std::size_t n_listed = 0;
      for (const auto& entry : entries) {
        n_listed += (bin_x(entry.first.max_corner().get<0>())
                     - bin_x(entry.first.min_corner().get<0>()) + 1)
            * (bin_y(entry.first.max_corner().get<1>())
               - bin_y(entry.first.min_corner().get<1>()) + 1);
      }
      if (n_listed <= 16 * entries.size() || (nx == 1 && ny == 1)) {
        break;
      }
      nx = std::max<std::size_t>(1, (nx + 1) / 2);
      ny = std::max<std::size_t>(1, (ny + 1) / 2);
    }

    // List the rectangles of each bin contiguously (compressed row storage)
    m_bin_begin.assign(m_nx * m_ny + 1, 0);
    for_each_bin_of(entries,
                    [&](std::size_t bin, const Entry&)
                    { ++m_bin_begin[bin + 1]; });
    for (std::size_t bin = 0; bin < m_nx * m_ny; ++bin) {
      m_bin_begin[bin + 1] += m_bin_begin[bin];
    }
    m_entries.resize(m_bin_begin.back());
    std::vector<std::size_t> next(m_bin_begin.begin(), m_bin_begin.end() - 1);
    for_each_bin_of(entries,
                    [&](std::size_t bin, const Entry& entry)
                    { m_entries[next[bin]++] = entry; });
  }

  /// Append the rectangles that contain the point to the candidates
  template<typename Container>
  void query(const Point& point, Container& candidates) const
  {
    const double x = point.get<0>();
    const double y = point.get<1>();
    // Outside the extent of all rectangles (or NaN): no rectangle contains it
    if (!(x >= m_x0 && x <= m_x1 && y >= m_y0 && y <= m_y1)) {
      return;
    }
    const std::size_t bin = bin_y(y) * m_nx + bin_x(x);
    for (std::size_t i = m_bin_begin[bin]; i < m_bin_begin[bin + 1]; ++i) {
      if (bg::intersects(m_entries[i].first, point)) {
        candidates.push_back(m_entries[i]);
      }
    }
  }

  /**
   * @brief Find the rectangles nearest to a point
   *
   * Searches rings of bins of growing size around the point's bin, until no
   * rectangle outside the searched bins can be nearer than the nearest found.
   * Appends all rectangles at the minimal distance (several only on an exact
   * tie) and returns true, or returns false without appending anything if the
   * point is not finite or the search would need more than max_rings rings.
   */
  template<typename Container>
  auto nearest(const Point& point, Container& candidates) const -> bool
  {
    // Gaps between neighboring cells are found within a ring or two. Points
    // in a large hole of a layer (such as between its two eta sides) are
    // found faster by the R-tree.
    constexpr std::size_t max_rings = 4;
    const double x = point.get<0>();
    const double y = point.get<1>();
    if (m_entries.empty() || !std::isfinite(x) || !std::isfinite(y)) {
      return false;
    }

    // Bin of the point, or of the nearest point of the grid's extent. The
    // comparisons come first, so that a distant point never converts a value
    // outside the range of std::size_t.
    const auto clamped_bin = [](double t, std::size_t n) -> std::size_t
    {
      if (!(t > 0)) {
        return 0;
      }
      if (t >= static_cast<double>(n - 1)) {
        return n - 1;
      }
      return static_cast<std::size_t>(t);
    };
    const std::size_t ix = clamped_bin((x - m_x0) * m_inv_wx, m_nx);
    const std::size_t iy = clamped_bin((y - m_y0) * m_inv_wy, m_ny);

    double best = std::numeric_limits<double>::infinity();
    const auto visit_bin = [&](std::size_t bx, std::size_t by)
    {
      const std::size_t bin = by * m_nx + bx;
      for (std::size_t i = m_bin_begin[bin]; i < m_bin_begin[bin + 1]; ++i) {
        const Entry& entry = m_entries[i];
        const double d = bg::comparable_distance(point, entry.first);
        if (d < best) {
          best = d;
          candidates.clear();
          candidates.push_back(entry);
        } else if (d == best
                   && std::none_of(candidates.begin(),
                                   candidates.end(),
                                   [&](const Entry& c)
                                   {
                                     return c.second == entry.second
                                         && bg::equals(c.first, entry.first);
                                   }))
        {
          // A rectangle is listed in every bin it overlaps: add it only once
          candidates.push_back(entry);
        }
      }
    };

    for (std::size_t r = 0; r <= max_rings; ++r) {
      const std::size_t ix_lo = ix >= r ? ix - r : 0;
      const std::size_t iy_lo = iy >= r ? iy - r : 0;
      const std::size_t ix_hi = std::min(ix + r, m_nx - 1);
      const std::size_t iy_hi = std::min(iy + r, m_ny - 1);

      // Visit the bins at ring distance r (the inner ones were visited)
      for (std::size_t by = iy_lo; by <= iy_hi; ++by) {
        if (by + r == iy || by == iy + r) {
          // Top or bottom row of the ring
          for (std::size_t bx = ix_lo; bx <= ix_hi; ++bx) {
            visit_bin(bx, by);
          }
        } else {
          // Left and right column of the ring, where inside the grid
          if (ix >= r) {
            visit_bin(ix - r, by);
          }
          if (ix + r < m_nx) {
            visit_bin(ix + r, by);
          }
        }
      }

      // Lower bound on the distance of the rectangles outside the visited
      // bins: they lie in the part of the extent beyond one of the sides of
      // the visited block. The sides are moved outwards by a margin that
      // covers the rounding of the bin computation.
      const bool all_visited =
          ix_lo == 0 && iy_lo == 0 && ix_hi == m_nx - 1 && iy_hi == m_ny - 1;
      if (all_visited) {
        return !candidates.empty();
      }
      double bound = std::numeric_limits<double>::infinity();
      const auto beyond =
          [&](double xmin, double ymin, double xmax, double ymax)
      {
        bound = std::min(bound,
                         bg::comparable_distance(
                             point, Box(Point(xmin, ymin), Point(xmax, ymax))));
      };
      if (ix_lo > 0) {
        beyond(m_x0, m_y0, edge_x(ix_lo) + m_margin_x, m_y1);
      }
      if (ix_hi < m_nx - 1) {
        beyond(edge_x(ix_hi + 1) - m_margin_x, m_y0, m_x1, m_y1);
      }
      if (iy_lo > 0) {
        beyond(m_x0, m_y0, m_x1, edge_y(iy_lo) + m_margin_y);
      }
      if (iy_hi < m_ny - 1) {
        beyond(m_x0, edge_y(iy_hi + 1) - m_margin_y, m_x1, m_y1);
      }
      if (bound > best) {
        return true;
      }
    }
    candidates.clear();
    return false;
  }

private:
  // Lower edge of a bin
  auto edge_x(std::size_t bx) const -> double
  {
    return m_x0 + static_cast<double>(bx) * m_wx;
  }
  auto edge_y(std::size_t by) const -> double
  {
    return m_y0 + static_cast<double>(by) * m_wy;
  }

  void set_bins(std::size_t nx, std::size_t ny)
  {
    m_nx = nx;
    m_ny = ny;
    m_inv_wx = m_x1 > m_x0 ? static_cast<double>(nx) / (m_x1 - m_x0) : 0;
    m_inv_wy = m_y1 > m_y0 ? static_cast<double>(ny) / (m_y1 - m_y0) : 0;
    m_wx = (m_x1 - m_x0) / static_cast<double>(nx);
    m_wy = (m_y1 - m_y0) / static_cast<double>(ny);
    // Far larger than the rounding of the bin computation, far smaller than
    // a bin
    const auto margin = [](double lo, double hi, double w)
    { return 1e-9 * std::max({std::abs(lo), std::abs(hi), w}); };
    m_margin_x = margin(m_x0, m_x1, m_wx);
    m_margin_y = margin(m_y0, m_y1, m_wy);
  }

  // Bin of a coordinate inside the extent; non-decreasing in the coordinate
  auto bin_x(double x) const -> std::size_t
  {
    return std::min(static_cast<std::size_t>((x - m_x0) * m_inv_wx), m_nx - 1);
  }
  auto bin_y(double y) const -> std::size_t
  {
    return std::min(static_cast<std::size_t>((y - m_y0) * m_inv_wy), m_ny - 1);
  }

  template<typename F>
  void for_each_bin_of(const std::vector<Entry>& entries, F&& f) const
  {
    for (const auto& entry : entries) {
      const std::size_t ix0 = bin_x(entry.first.min_corner().get<0>());
      const std::size_t ix1 = bin_x(entry.first.max_corner().get<0>());
      const std::size_t iy0 = bin_y(entry.first.min_corner().get<1>());
      const std::size_t iy1 = bin_y(entry.first.max_corner().get<1>());
      for (std::size_t iy = iy0; iy <= iy1; ++iy) {
        for (std::size_t ix = ix0; ix <= ix1; ++ix) {
          f(iy * m_nx + ix, entry);
        }
      }
    }
  }

  double m_x0 {0}, m_y0 {0}, m_x1 {-1}, m_y1 {-1};
  double m_inv_wx {0}, m_inv_wy {0};
  double m_wx {0}, m_wy {0};
  double m_margin_x {0}, m_margin_y {0};
  std::size_t m_nx {1}, m_ny {1};
  std::vector<std::size_t> m_bin_begin;
  std::vector<Entry> m_entries;
};
}  // namespace

struct RTree::Impl
{
  Tree tree;
  Grid grid;
  // Cell ID of each position in the constructor's cells
  std::vector<uint64_t> ids;
};

RTree::RTree(RTreeHelpers::CoordinateSystem coord_sys,
             const std::vector<Cell>& cells)
    : m_coordinate_system(coord_sys)
{
  const bool is_xyz = coord_sys == RTreeHelpers::CoordinateSystem::XYZ;

  std::vector<Entry> entries;
  std::vector<uint64_t> ids;
  ids.reserve(cells.size());
  for (const Cell& cell : cells) {
    // One rectangle per cell, or two if the cell crosses phi = +-pi
    const auto rectangles = is_xyz
        ? std::vector {RTreeHelpers::build_xy_box(
              cell.x(), cell.y(), cell.dx(), cell.dy())}
        : RTreeHelpers::build_eta_phi_boxes(
              cell.eta(), cell.phi(), cell.deta(), cell.dphi());

    for (const auto& r : rectangles) {
      entries.emplace_back(Box(Point(r[0], r[1]), Point(r[2], r[3])),
                           ids.size());
    }
    ids.push_back(cell.id());
  }

  // Constructing the tree from all rectangles at once packs it (bulk loading)
  m_impl = std::make_unique<const Impl>(Impl {
      Tree(entries.begin(), entries.end()), Grid(entries), std::move(ids)});
}

RTree::~RTree() = default;

auto RTree::query_point(const Position& pos) const -> uint64_t
{
  return m_impl->ids[query_point_index(pos)];
}

auto RTree::query_point_index(const Position& pos) const -> std::size_t
{
  const Point point = m_coordinate_system == RTreeHelpers::CoordinateSystem::XYZ
      ? Point(pos.x(), pos.y())
      : Point(pos.eta(), pos.phi());

  // The cells whose rectangle contains the point: normally exactly one, but
  // the rectangles of neighboring cells can overlap slightly
  boost::container::small_vector<Entry, 4> candidates;
  m_impl->grid.query(point, candidates);

  // In a gap between cells or outside the layer: the nearest rectangle. The
  // grid finds it unless the point is far from all rectangles.
  if (candidates.empty() && !m_impl->grid.nearest(point, candidates)) {
    // The R-tree yields the rectangles in order of distance: take all at the
    // minimal distance, so that the rule below also decides ties here
    const Tree& tree = m_impl->tree;
    double min_distance = 0;
    for (auto it = tree.qbegin(bgi::nearest(point, tree.size()));
         it != tree.qend();
         ++it)
    {
      const double d = bg::comparable_distance(point, it->first);
      if (!candidates.empty() && d != min_distance) {
        break;
      }
      min_distance = d;
      candidates.push_back(*it);
      // A point that is not finite has no meaningful nearest rectangle
      if (!std::isfinite(d)) {
        break;
      }
    }
  }
  if (candidates.empty()) {
    throw std::invalid_argument("No cell found for query point");
  }

  // Of several candidates, take the cell whose rectangle center is nearest
  // (and the smaller cell ID on an exact tie)
  const auto& ids = m_impl->ids;
  const auto center_distance_and_id = [&point, &ids](const Entry& entry)
  {
    const auto center = bg::return_centroid<Point>(entry.first);
    return std::pair(bg::comparable_distance(point, center), ids[entry.second]);
  };
  const auto best = std::min_element(
      candidates.begin(),
      candidates.end(),
      [&](const Entry& a, const Entry& b)
      { return center_distance_and_id(a) < center_distance_and_id(b); });
  return best->second;
}
