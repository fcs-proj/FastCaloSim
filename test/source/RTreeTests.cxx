// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

// Tests of the cell lookup (RTree): first the lookup rule on small hand-made
// layers, then a comparison with a brute-force search on the ATLAS geometry.

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <random>
#include <string>
#include <utility>
#include <vector>

#include "FastCaloSim/Geometry/RTree.h"

#include <gtest/gtest.h>

#include "FastCaloSim/Geometry/CaloGeo.h"
#include "ROOT/RDataFrame.hxx"
#include "TestConfig/AtlasGeoTestsConfig.h"

namespace
{
const auto kEtaPhi = RTreeHelpers::CoordinateSystem::EtaPhiR;

// A cell centred at (eta, phi) with full widths (deta, dphi)
auto eta_phi_cell(uint64_t id, double eta, double phi, double deta, double dphi)
    -> Cell
{
  const Position centre {0, 0, 0, eta, phi, 0};
  return Cell(
      id, centre, 0, true, false, true, false, false, 0, 0, 0, deta, dphi, 0);
}

auto at_eta_phi(double eta, double phi) -> Position
{
  return Position {0, 0, 0, eta, phi, 0};
}
}  // namespace

TEST(RTreeTests, PositionInsideOneCell)
{
  const RTree rtree(kEtaPhi,
                    {eta_phi_cell(1, 0.5, 0.0, 1.0, 0.2),  // eta in [0, 1]
                     eta_phi_cell(2, 1.5, 0.0, 1.0, 0.2)});  // eta in [1, 2]

  EXPECT_EQ(rtree.query_point(at_eta_phi(0.3, 0.05)), 1);
  EXPECT_EQ(rtree.query_point(at_eta_phi(1.7, -0.05)), 2);
}

TEST(RTreeTests, PositionInOverlapOfTwoCells)
{
  // The cells overlap for eta in [0.875, 1.125]
  const RTree rtree(
      kEtaPhi,
      {eta_phi_cell(1, 0.5, 0.0, 1.25, 0.2),  // eta in [-0.125, 1.125]
       eta_phi_cell(2, 1.5, 0.0, 1.25, 0.2)});  // eta in [0.875, 2.125]

  // The cell with the nearest centre is taken
  EXPECT_EQ(rtree.query_point(at_eta_phi(0.9, 0.0)), 1);
  EXPECT_EQ(rtree.query_point(at_eta_phi(1.1, 0.0)), 2);
  // Both centres equally far: the smaller ID
  EXPECT_EQ(rtree.query_point(at_eta_phi(1.0, 0.0)), 1);
}

TEST(RTreeTests, PositionOutsideAllCells)
{
  // There is a gap for eta in [1, 2]
  const RTree rtree(kEtaPhi,
                    {eta_phi_cell(1, 0.5, 0.0, 1.0, 0.2),  // eta in [0, 1]
                     eta_phi_cell(2, 2.5, 0.0, 1.0, 0.2)});  // eta in [2, 3]

  // The cell with the nearest rectangle is taken
  EXPECT_EQ(rtree.query_point(at_eta_phi(1.2, 0.0)), 1);
  EXPECT_EQ(rtree.query_point(at_eta_phi(1.8, 0.0)), 2);
  EXPECT_EQ(rtree.query_point(at_eta_phi(-1000, 0.0)), 1);
  EXPECT_EQ(rtree.query_point(at_eta_phi(1000, 0.0)), 2);
}

TEST(RTreeTests, CellAcrossPhiBoundary)
{
  // Cell 1 covers phi in [pi - 0.1, pi] and [-pi, -pi + 0.1]
  const RTree rtree(kEtaPhi,
                    {eta_phi_cell(1, 0.5, M_PI, 1.0, 0.2),
                     eta_phi_cell(2, 0.5, 0.0, 1.0, 0.2)});

  EXPECT_EQ(rtree.query_point(at_eta_phi(0.5, M_PI - 0.05)), 1);
  EXPECT_EQ(rtree.query_point(at_eta_phi(0.5, -M_PI + 0.05)), 1);
  EXPECT_EQ(rtree.query_point(at_eta_phi(0.5, 0.05)), 2);
}

// ---- Comparison with a brute-force search on the ATLAS geometry

namespace
{
// A cell rectangle in lookup coordinates (u, v): (x, y) or (eta, phi)
struct Rect
{
  double u_min, v_min, u_max, v_max;
  uint64_t cell_id;

  auto contains(double u, double v) const -> bool
  {
    return u_min <= u && u <= u_max && v_min <= v && v <= v_max;
  }

  auto distance2(double u, double v) const -> double
  {
    const double du = std::max({u_min - u, 0., u - u_max});
    const double dv = std::max({v_min - v, 0., v - v_max});
    return du * du + dv * dv;
  }

  auto centre_distance2(double u, double v) const -> double
  {
    const double du = 0.5 * (u_min + u_max) - u;
    const double dv = 0.5 * (v_min + v_max) - v;
    return du * du + dv * dv;
  }
};
}  // namespace

class AtlasRTreeTests : public ::testing::Test
{
protected:
  inline static CaloGeo* geo = nullptr;

  static void SetUpTestSuite()
  {
    ROOT::RDataFrame df(AtlasGeoTestsConfig::ATLAS_CALO_CELL_TREE_NAME,
                        AtlasGeoTestsConfig::ATLAS_CALO_CELL_PATH);
    char tmpl[] = "/tmp/fastcalosimXXXXXX";
    const std::string geo_dir = mkdtemp(tmpl);
    geo = new CaloGeo();
    geo->build(df, geo_dir);
    geo->load(geo_dir);
  }

  // The rectangles of all cells of a layer
  static auto cell_rects(unsigned int layer) -> std::vector<Rect>
  {
    std::vector<Rect> rects;
    for (unsigned int i = 0; i < geo->n_cells(layer); ++i) {
      const Cell cell = geo->get_cell_at_idx(layer, i);
      const auto boxes = geo->is_xyz(layer)
          ? std::vector {RTreeHelpers::build_xy_box(
                cell.x(), cell.y(), cell.dx(), cell.dy())}
          : RTreeHelpers::build_eta_phi_boxes(
                cell.eta(), cell.phi(), cell.deta(), cell.dphi());
      for (const auto& b : boxes) {
        rects.push_back({b[0], b[1], b[2], b[3], cell.id()});
      }
    }
    return rects;
  }

  // Is this the cell the lookup rule gives for (u, v)? Checked against all
  // rectangles of the layer.
  static auto is_expected_cell(uint64_t cell_id,
                               const std::vector<Rect>& rects,
                               double u,
                               double v) -> bool
  {
    // Of the rectangles containing (u, v): the nearest centre, then smaller ID
    const Rect* best = nullptr;
    for (const auto& r : rects) {
      if (!r.contains(u, v)) {
        continue;
      }
      if (best == nullptr
          || std::pair(r.centre_distance2(u, v), r.cell_id)
              < std::pair(best->centre_distance2(u, v), best->cell_id))
      {
        best = &r;
      }
    }
    if (best != nullptr) {
      return cell_id == best->cell_id;
    }

    // No rectangle contains (u, v): a cell with a rectangle nearest to it
    double min_distance2 = std::numeric_limits<double>::infinity();
    for (const auto& r : rects) {
      min_distance2 = std::min(min_distance2, r.distance2(u, v));
    }
    for (const auto& r : rects) {
      if (r.cell_id == cell_id && r.distance2(u, v) == min_distance2) {
        return true;
      }
    }
    return false;
  }
};

TEST_F(AtlasRTreeTests, SameCellAsBruteForceSearch)
{
  std::mt19937_64 rng(42);
  std::normal_distribution<double> gauss(0., 1.);
  const double inf = std::numeric_limits<double>::infinity();

  for (unsigned int layer = 0; layer < geo->n_layers(); ++layer) {
    const auto rects = cell_rects(layer);

    // Positions around 300 random cells: scattered inside and around the
    // cell, far outside the layer, on the cell's corners and edges, and one
    // representable number on either side of each edge
    std::vector<std::array<double, 2>> points;
    std::uniform_int_distribution<size_t> pick(0, rects.size() - 1);
    for (int i = 0; i < 300; ++i) {
      const Rect& r = rects[pick(rng)];
      const double u_mid = 0.5 * (r.u_min + r.u_max);
      const double v_mid = 0.5 * (r.v_min + r.v_max);
      const double du = r.u_max - r.u_min;
      const double dv = r.v_max - r.v_min;

      for (int k = 0; k < 10; ++k) {
        points.push_back(
            {u_mid + 2 * du * gauss(rng), v_mid + 2 * dv * gauss(rng)});
      }
      points.push_back({u_mid + 50 * du, v_mid});
      for (double u : {r.u_min, u_mid, r.u_max}) {
        for (double v : {r.v_min, v_mid, r.v_max}) {
          points.push_back({u, v});
        }
      }
      for (double u : {r.u_min, r.u_max}) {
        points.push_back({std::nextafter(u, -inf), v_mid});
        points.push_back({std::nextafter(u, inf), v_mid});
      }
      for (double v : {r.v_min, r.v_max}) {
        points.push_back({u_mid, std::nextafter(v, -inf)});
        points.push_back({u_mid, std::nextafter(v, inf)});
      }
    }

    for (const auto& [u, v] : points) {
      const Position pos = geo->is_xyz(layer) ? Position {u, v, 0, 0, 0, 0}
                                              : Position {0, 0, 0, u, v, 0};
      const uint64_t cell_id = geo->get_cell(layer, pos).id();

      ASSERT_TRUE(is_expected_cell(cell_id, rects, u, v))
          << "layer " << layer << ", position (" << u << ", " << v << "), cell "
          << cell_id;
    }
  }
}
