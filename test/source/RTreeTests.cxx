// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

// Tests of the cell lookup (RTree): first the lookup rule on small hand-made
// layers, then a comparison with a brute-force search on the ATLAS geometry.

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <functional>
#include <limits>
#include <memory>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "FastCaloSim/Geometry/RTree.h"

#include <gtest/gtest.h>

#include "FastCaloSim/Core/TFCSParametrizationBase.h"
#include "FastCaloSim/Core/TFCSSimulationState.h"
#include "FastCaloSim/Definitions/FloatCompare.h"
#include "FastCaloSim/Definitions/Phi.h"
#include "FastCaloSim/Geometry/CaloGeo.h"
#include "ROOT/RDataFrame.hxx"
#include "TClass.h"
#include "TTree.h"
#include "TestConfig/AtlasGeoTestsConfig.h"

namespace
{
const auto kEtaPhi = RTreeHelpers::CoordinateSystem::EtaPhiR;

// A cell centered at (eta, phi) with full widths (deta, dphi)
auto eta_phi_cell(uint64_t id, double eta, double phi, double deta, double dphi)
    -> Cell
{
  const Position center {
      .m_x = 0, .m_y = 0, .m_z = 0, .m_eta = eta, .m_phi = phi, .m_r = 0};
  return {id,
          center,
          0,
          /*isBarrel=*/true,
          /*isXYZ=*/false,
          /*isEtaPhiR=*/true,
          /*isEtaPhiZ=*/false,
          /*isRPhiZ=*/false,
          0,
          0,
          0,
          deta,
          dphi,
          0};
}

auto at_eta_phi(double eta, double phi) -> Position
{
  return Position {
      .m_x = 0, .m_y = 0, .m_z = 0, .m_eta = eta, .m_phi = phi, .m_r = 0};
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

  // The cell with the nearest center is taken
  EXPECT_EQ(rtree.query_point(at_eta_phi(0.9, 0.0)), 1);
  EXPECT_EQ(rtree.query_point(at_eta_phi(1.1, 0.0)), 2);
  // Both centers equally far: the smaller ID
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

TEST(RTreeTests, PositionNotFinite)
{
  const RTree rtree(kEtaPhi,
                    {eta_phi_cell(1, 0.5, 0.0, 1.0, 0.2),  // eta in [0, 1]
                     eta_phi_cell(2, 1.5, 0.0, 1.0, 0.2)});  // eta in [1, 2]

  // A position that is not finite lies in no cell; one of the cells is
  // returned, without undefined behavior in the cell lookup
  const double nan = std::numeric_limits<double>::quiet_NaN();
  const double inf = std::numeric_limits<double>::infinity();
  for (const Position& pos : {at_eta_phi(nan, 0.0),
                              at_eta_phi(0.5, nan),
                              at_eta_phi(nan, nan),
                              at_eta_phi(inf, 0.0),
                              at_eta_phi(0.5, -inf)})
  {
    const auto id = rtree.query_point(pos);
    EXPECT_TRUE(id == 1 || id == 2) << "id " << id;
  }
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

TEST(RTreeTests, AzimuthOutsideMinusPiToPi)
{
  // A hit placed as center + offset can have an azimuth just beyond +-pi: it
  // must be found at the same direction inside [-pi, pi], not in the nearest
  // cell on the other side of the seam
  const RTree rtree(kEtaPhi,
                    {eta_phi_cell(1, 0.5, M_PI - 0.05, 1.0, 0.1),
                     eta_phi_cell(2, 0.5, -M_PI + 0.05, 1.0, 0.1)});

  EXPECT_EQ(rtree.query_point(at_eta_phi(0.5, M_PI + 0.03)), 2);
  EXPECT_EQ(rtree.query_point(at_eta_phi(0.5, -M_PI - 0.03)), 1);
  EXPECT_EQ(rtree.query_point(at_eta_phi(0.5, M_PI - 0.03)), 1);
  EXPECT_EQ(rtree.query_point(at_eta_phi(0.5, -M_PI + 0.03)), 2);
}

TEST(RTreeTests, PositionFarOutsideAllCells)
{
  const RTree rtree(kEtaPhi,
                    {eta_phi_cell(1, 0.5, 0.0, 1.0, 0.2),  // eta in [0, 1]
                     eta_phi_cell(2, 2.5, 0.0, 1.0, 0.2)});  // eta in [2, 3]

  // So far away that the distances to both cells are equal in floating
  // point, and so are those to their centers: the smaller ID
  EXPECT_EQ(rtree.query_point(at_eta_phi(-1e30, 0.0)), 1);
  EXPECT_EQ(rtree.query_point(at_eta_phi(1e30, 0.0)), 1);
  EXPECT_EQ(rtree.query_point(at_eta_phi(1e30, 1e30)), 1);
}

TEST(RTreeTests, TieFarFromAllCells)
{
  // Unit cells at eta = 0..79 and 221..300. The position at eta = 150 is far
  // from both groups and equally near to the cells at eta = 79 and 221, whose
  // centers are also equally far: the smaller ID is taken, whichever side
  for (const bool left_smaller : {true, false}) {
    std::vector<Cell> cells;
    for (int i = 0; i <= 300; ++i) {
      if (i >= 80 && i <= 220) {
        continue;
      }
      const auto id = static_cast<uint64_t>(
          (i < 150) == left_smaller ? 1000 + i : 2000 + i);
      cells.push_back(eta_phi_cell(id, i, 0.0, 1.0, 0.2));
    }
    const RTree rtree(kEtaPhi, cells);
    EXPECT_EQ(rtree.query_point(at_eta_phi(150, 0.0)),
              left_smaller ? 1079 : 1221);
  }
}

TEST(RTreeTests, ManyNearestCellsTied)
{
  // Cells 7, 3, 9 and 5 have the same rectangle, at eta in [0, 1]; the other
  // cells are at eta >= 200. The position at eta = 60 is equally near to the
  // four (and so are their centers): the smallest ID
  std::vector<Cell> cells;
  for (const uint64_t id : std::array<uint64_t, 4> {7, 3, 9, 5}) {
    cells.push_back(eta_phi_cell(id, 0.5, 0.0, 1.0, 0.2));
  }
  for (uint64_t i = 200; i <= 300; ++i) {
    cells.push_back(
        eta_phi_cell(1000 + i, static_cast<double>(i), 0.0, 1.0, 0.2));
  }
  const RTree rtree(kEtaPhi, cells);
  EXPECT_EQ(rtree.query_point(at_eta_phi(60, 0.0)), 3);
  EXPECT_EQ(rtree.query_point(at_eta_phi(60, 0.05)), 3);
  EXPECT_EQ(rtree.query_point(at_eta_phi(-60, 1.0)), 3);
}

TEST(RTreeTests, CellsOfVeryDifferentWidths)
{
  // The median cell width is far smaller than the extent of the layer
  const RTree rtree(kEtaPhi,
                    {eta_phi_cell(1, 0.0, 0.0, 1e-20, 0.2),
                     eta_phi_cell(2, 1e-19, 0.0, 1e-20, 0.2),
                     eta_phi_cell(3, 1.0, 0.0, 0.1, 0.2)});

  EXPECT_EQ(rtree.query_point(at_eta_phi(1e-19, 0.0)), 2);
  EXPECT_EQ(rtree.query_point(at_eta_phi(1.0, 0.0)), 3);
  EXPECT_EQ(rtree.query_point(at_eta_phi(0.5, 0.0)), 3);
}

TEST(RTreeTests, LargeLayerIds)
{
  // Layer IDs are not limited to a small, dense range
  class FixedCellGeo : public CaloGeo
  {
  public:
    explicit FixedCellGeo(uint64_t id)
        : m_cell(eta_phi_cell(id, 0.0, 0.0, 1.0, 1.0))
    {
    }
    auto get_cell(std::size_t /*layer*/, const Position& /*pos*/) const
        -> Cell override
    {
      return m_cell;
    }

  private:
    Cell m_cell;
  };

  CaloGeo geo;
  const std::size_t max_layer = std::numeric_limits<std::size_t>::max();
  geo.set_alt_geo_handler(7, std::make_shared<FixedCellGeo>(1));
  geo.set_alt_geo_handler(100000000, std::make_shared<FixedCellGeo>(2));
  geo.set_alt_geo_handler(max_layer, std::make_shared<FixedCellGeo>(3));

  EXPECT_EQ(geo.get_cell(7, at_eta_phi(0, 0)).id(), 1);
  EXPECT_EQ(geo.get_cell(100000000, at_eta_phi(0, 0)).id(), 2);
  EXPECT_EQ(geo.get_cell(max_layer, at_eta_phi(0, 0)).id(), 3);
  EXPECT_THROW(geo.get_cell(0, at_eta_phi(0, 0)), std::runtime_error);
  EXPECT_THROW(geo.get_cell(8, at_eta_phi(0, 0)), std::runtime_error);
  EXPECT_THROW(geo.get_cell(100000001, at_eta_phi(0, 0)), std::runtime_error);
}

namespace
{
// A cell of a test geometry: id, coordinate system (one of "XYZ", "EtaPhiR",
// "EtaPhiZ", "RPhiZ") and layer
struct TestCell
{
  unsigned long long id;
  std::string system;
  unsigned int layer;
};

// Builds a geometry from the cells, then loads it and runs check on it
void build_geometry(const std::vector<TestCell>& cells,
                    const std::function<void(CaloGeo&)>& check = {})
{
  unsigned int layer = 0;
  bool is_barrel = true;
  unsigned long long id = 0;
  float x = 1000.F;
  float y = 0.F;
  float z = 0.F;
  float phi = 0.F;
  float eta = 0.F;
  float r = 1000.F;
  float dx = 0.F;
  float dy = 0.F;
  float dz = 10.F;
  float dphi = 0.1F;
  float deta = 0.1F;
  float dr = 10.F;
  bool is_xyz = false;
  bool is_eta_phi_r = false;
  bool is_eta_phi_z = false;
  bool is_r_phi_z = false;
  TTree tree("cells", "cells");
  tree.SetDirectory(nullptr);
  tree.Branch("layer", &layer);
  tree.Branch("isBarrel", &is_barrel);
  tree.Branch("id", &id);
  tree.Branch("x", &x);
  tree.Branch("y", &y);
  tree.Branch("z", &z);
  tree.Branch("phi", &phi);
  tree.Branch("eta", &eta);
  tree.Branch("r", &r);
  tree.Branch("dx", &dx);
  tree.Branch("dy", &dy);
  tree.Branch("dz", &dz);
  tree.Branch("dphi", &dphi);
  tree.Branch("deta", &deta);
  tree.Branch("dr", &dr);
  tree.Branch("isXYZ", &is_xyz);
  tree.Branch("isEtaPhiR", &is_eta_phi_r);
  tree.Branch("isEtaPhiZ", &is_eta_phi_z);
  tree.Branch("isRPhiZ", &is_r_phi_z);
  for (const auto& [cell_id, system, cell_layer] : cells) {
    id = cell_id;
    layer = cell_layer;
    eta = 0.2F * static_cast<float>(cell_id);
    is_xyz = system == "XYZ";
    is_eta_phi_r = system == "EtaPhiR";
    is_eta_phi_z = system == "EtaPhiZ";
    is_r_phi_z = system == "RPhiZ";
    tree.Fill();
  }
  ROOT::RDataFrame geo_df(tree);
  // build() writes the cell store into this directory, one per test, as ctest
  // runs the tests in parallel
  const std::filesystem::path dir = std::filesystem::temp_directory_path()
      / (std::string("fcs_rtree_tests_geo_")
         + ::testing::UnitTest::GetInstance()->current_test_info()->name());
  std::filesystem::create_directories(dir);
  CaloGeo geo;
  try {
    geo.build(geo_df, dir.string());
    if (check) {
      // The R-trees for the cell lookup are built when loading
      geo.load(dir.string());
      check(geo);
    }
  } catch (...) {
    std::filesystem::remove_all(dir);
    throw;
  }
  std::filesystem::remove_all(dir);
}
}  // namespace

TEST(RTreeTests, CellsInRPhiZAreRejected)
{
  // The cell lookup works in (eta, phi) or (x, y); a geometry with cells
  // described in (r, phi, z) must be rejected instead of being looked up with
  // an eta extent it does not have, also when the cell is not the first of
  // its layer
  EXPECT_NO_THROW(build_geometry({{1, "EtaPhiR", 3}, {2, "EtaPhiR", 3}}));
  EXPECT_THROW(build_geometry({{1, "RPhiZ", 3}}), std::runtime_error);
  EXPECT_THROW(build_geometry({{1, "EtaPhiR", 3}, {2, "RPhiZ", 3}}),
               std::runtime_error);
}

TEST(RTreeTests, SparseLayerIds)
{
  // Layer IDs need not start at 0 or be contiguous: the geometry lists the
  // IDs of its layers, which the simulation loops over
  build_geometry(
      {{.id = 1, .system = "EtaPhiR", .layer = 3},
       {.id = 2, .system = "EtaPhiR", .layer = 3},
       {.id = 3, .system = "EtaPhiR", .layer = 7}},
      [](CaloGeo& geo) -> void
      {
        EXPECT_EQ(geo.n_layers(), 2U);
        EXPECT_EQ(geo.layer_ids(), (std::vector<std::size_t> {3, 7}));
        EXPECT_EQ(geo.layer_index(3), 0U);
        EXPECT_EQ(geo.layer_index(7), 1U);
        EXPECT_THROW(static_cast<void>(geo.layer_index(5)), std::out_of_range);
        EXPECT_TRUE(geo.has_layer(3));
        EXPECT_TRUE(geo.has_layer(7));
        EXPECT_FALSE(geo.has_layer(0));
        EXPECT_FALSE(geo.has_layer(5));
        EXPECT_FALSE(geo.has_layer(8));
        EXPECT_EQ(geo.get_cell(7, at_eta_phi(0.6, 0)).id(), 3U);
      });
}

TEST(RTreeTests, EnergyRenormalizationWithDistantLayerIds)
{
  // Per-layer data has one entry per layer, also for layer IDs far apart
  constexpr std::size_t far_layer = 100000000;
  build_geometry(
      {{.id = 1, .system = "EtaPhiR", .layer = 3},
       {.id = 2, .system = "EtaPhiR", .layer = far_layer}},
      [](CaloGeo& geo) -> void
      {
        EXPECT_EQ(geo.layer_ids(), (std::vector<std::size_t> {3, far_layer}));
        EXPECT_EQ(geo.layer_index(far_layer), 1U);
        // The class is not exported: create it through ROOT and own it
        // through the exported base class
        const std::unique_ptr<TFCSParametrizationBase> renormalization(
            static_cast<TFCSParametrizationBase*>(
                TClass::GetClass("TFCSEnergyRenormalization")->New()));
        renormalization->set_geometry(&geo);
        TFCSSimulationState state;
        state.set_E(3, 10.);
        state.set_E(far_layer, 20.);
        state.deposit(1, 5.);
        state.deposit(2, 40.);
        EXPECT_EQ(renormalization->simulate(state, nullptr, nullptr),
                  FCSSuccess);
        // The cells are scaled to the energies of their layers
        EXPECT_DOUBLE_EQ(state.cells().at(1), 10.);
        EXPECT_DOUBLE_EQ(state.cells().at(2), 20.);
      });
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

  auto center_distance2(double u, double v) const -> double
  {
    const double du = 0.5 * (u_min + u_max) - u;
    const double dv = 0.5 * (v_min + v_max) - v;
    return du * du + dv * dv;
  }
};

class AtlasRTreeTests : public ::testing::Test
{
protected:
  inline static CaloGeo* geo = nullptr;

  static void SetUpTestSuite()
  {
    ROOT::RDataFrame df(AtlasGeoTestsConfig::ATLAS_CALO_CELL_TREE_NAME,
                        AtlasGeoTestsConfig::ATLAS_CALO_CELL_PATH);
    std::string tmpl = "/tmp/fastcalosimXXXXXX";
    const char* dir = mkdtemp(tmpl.data());
    if (dir == nullptr) {
      throw std::runtime_error("Cannot create a temporary directory");
    }
    const std::string geo_dir = dir;
    geo = new CaloGeo();
    geo->build(df, geo_dir);
    geo->load(geo_dir);
  }

  // The rectangles of all cells of a layer
  static auto cell_rects(std::size_t layer) -> std::vector<Rect>
  {
    std::vector<Rect> rects;
    for (std::size_t i = 0; i < geo->n_cells(layer); ++i) {
      const Cell cell = geo->get_cell_at_idx(layer, i);
      const auto boxes = geo->is_xyz(layer)
          ? std::vector {RTreeHelpers::build_xy_box(
                cell.x(), cell.y(), cell.dx(), cell.dy())}
          : RTreeHelpers::build_eta_phi_boxes(
                cell.eta(), cell.phi(), cell.deta(), cell.dphi());
      for (const auto& b : boxes) {
        rects.push_back({.u_min = b[0],
                         .v_min = b[1],
                         .u_max = b[2],
                         .v_max = b[3],
                         .cell_id = cell.id()});
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
    // Of the rectangles containing (u, v): the nearest center, then smaller ID
    const Rect* best = nullptr;
    for (const auto& r : rects) {
      if (!r.contains(u, v)) {
        continue;
      }
      if (best == nullptr
          || std::pair(r.center_distance2(u, v), r.cell_id)
              < std::pair(best->center_distance2(u, v), best->cell_id))
      {
        best = &r;
      }
    }
    if (best != nullptr) {
      return cell_id == best->cell_id;
    }

    // No rectangle contains (u, v): of the rectangles nearest to it, the
    // nearest center, then the smaller ID
    double min_distance2 = std::numeric_limits<double>::infinity();
    for (const auto& r : rects) {
      min_distance2 = std::min(min_distance2, r.distance2(u, v));
    }
    for (const auto& r : rects) {
      if (!FastCaloSim::exactly_equal(r.distance2(u, v), min_distance2)) {
        continue;
      }
      if (best == nullptr
          || std::pair(r.center_distance2(u, v), r.cell_id)
              < std::pair(best->center_distance2(u, v), best->cell_id))
      {
        best = &r;
      }
    }
    return best != nullptr && cell_id == best->cell_id;
  }
};
}  // namespace

TEST_F(AtlasRTreeTests, SameCellAsBruteForceSearch)
{
  // A fixed seed, for the same test positions in every run
  // NOLINTBEGIN(bugprone-random-generator-seed,cert-msc32-c,cert-msc51-cpp)
  std::mt19937_64 rng(42);
  // NOLINTEND(bugprone-random-generator-seed,cert-msc32-c,cert-msc51-cpp)
  std::normal_distribution<double> gauss(0., 1.);
  const double inf = std::numeric_limits<double>::infinity();

  for (const std::size_t layer : geo->layer_ids()) {
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
      for (double const u : {r.u_min, u_mid, r.u_max}) {
        for (double const v : {r.v_min, v_mid, r.v_max}) {
          points.push_back({u, v});
        }
      }
      for (double const u : {r.u_min, r.u_max}) {
        points.push_back({std::nextafter(u, -inf), v_mid});
        points.push_back({std::nextafter(u, inf), v_mid});
      }
      for (double const v : {r.v_min, r.v_max}) {
        points.push_back({u_mid, std::nextafter(v, -inf)});
        points.push_back({u_mid, std::nextafter(v, inf)});
      }
    }

    // Positions anywhere around the layer, many of them in holes of the
    // layer or beyond its edges, where the R-tree finds the nearest cell
    double u_lo = std::numeric_limits<double>::infinity();
    double v_lo = u_lo;
    double u_hi = -u_lo;
    double v_hi = -u_lo;
    for (const auto& r : rects) {
      u_lo = std::min(u_lo, r.u_min);
      v_lo = std::min(v_lo, r.v_min);
      u_hi = std::max(u_hi, r.u_max);
      v_hi = std::max(v_hi, r.v_max);
    }
    std::uniform_real_distribution<double> u_around(u_lo - 0.5 * (u_hi - u_lo),
                                                    u_hi + 0.5 * (u_hi - u_lo));
    std::uniform_real_distribution<double> v_around(v_lo - 0.5 * (v_hi - v_lo),
                                                    v_hi + 0.5 * (v_hi - v_lo));
    for (int i = 0; i < 1000; ++i) {
      points.push_back({u_around(rng), v_around(rng)});
    }

    for (const auto& [u, v] : points) {
      const Position pos = geo->is_xyz(layer)
          ? Position {.m_x = u,
                      .m_y = v,
                      .m_z = 0,
                      .m_eta = 0,
                      .m_phi = 0,
                      .m_r = 0}
          : Position {
                .m_x = 0, .m_y = 0, .m_z = 0, .m_eta = u, .m_phi = v, .m_r = 0};
      const uint64_t cell_id = geo->get_cell(layer, pos).id();

      // An azimuth outside [-pi, pi] is the same direction as its wrapped
      // value, which is where the cell rectangles are
      const double v_ref = geo->is_xyz(layer) ? v : FastCaloSim::wrap_phi(v);
      ASSERT_TRUE(is_expected_cell(cell_id, rects, u, v_ref))
          << "layer " << layer << ", position (" << u << ", " << v << "), cell "
          << cell_id;
    }
  }
}
