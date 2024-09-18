#include "AtlasGeoTests.h"

#include <gtest/gtest.h>

#include "FastCaloSim/Geometry/CaloGeo.h"

TEST_F(AtlasGeoTests, GeometryHealth)
{
  // Check that the barrel layers are indeed barrel layers
  std::vector<int> barrel_layers = {0, 1, 2, 3, 12, 13, 14, 15, 16, 18, 19, 20};
  for (int layer : barrel_layers) {
    ASSERT_TRUE(AtlasGeoTests::geo->is_barrel(layer));
  }
  // Check that the endcap layers are indeed endcap layers
  std::vector<int> endcap_layers = {4, 5, 6, 7, 8, 9, 10, 11, 17, 21, 22, 23};
  for (int layer : endcap_layers) {
    ASSERT_FALSE(AtlasGeoTests::geo->is_barrel(layer));
  }
}

struct Hit2D
{
  auto x() const -> double { return m_x; }
  auto y() const -> double { return m_y; }
  auto eta() const -> double { return m_eta; }
  auto phi() const -> double { return m_phi; }
  double m_x, m_y, m_z, m_eta, m_phi, m_r;
};

TEST_F(AtlasGeoTests, PreSamplerBLookup)
{
  constexpr int layer = 0;

  // Loop over all cells in the layer
  for (int i = 0; i < AtlasGeoTests::geo->n_cells(layer); ++i) {
    const auto& cell = AtlasGeoTests::geo->get_cell_at_idx(layer, i);

    // Loop over random points within the cell
    for (int irnd = 0; irnd < AtlasGeoTestsConfig::N_RANDOM_POINTS; ++irnd) {
      Hit2D hit {};
      hit.m_eta = AtlasGeoTestsConfig::sample(cell.eta() - 0.49 * cell.deta(),
                                              cell.eta() + 0.49 * cell.deta());
      hit.m_phi = Cell::norm_angle(AtlasGeoTestsConfig::sample(
          cell.phi() - 0.49 * cell.dphi(), cell.phi() + 0.49 * cell.dphi()));

      // Retrieve the cell closest to the hit position
      const auto& best_cell = AtlasGeoTests::geo->get_cell(layer, hit);

      // Check if the hit position is inside the cell
      ASSERT_NEAR(hit.m_eta, best_cell.eta(), 0.5 * best_cell.deta());
      ASSERT_NEAR(Cell::norm_angle(hit.m_phi - best_cell.phi()),
                  0.0,
                  0.5 * best_cell.dphi());
    }
  }
}

TEST_F(AtlasGeoTests, EME1Lookup)
{
  constexpr int layer = 5;

  // Loop over all cells in the layer
  for (int i = 0; i < AtlasGeoTests::geo->n_cells(layer); ++i) {
    const auto& cell = AtlasGeoTests::geo->get_cell_at_idx(layer, i);

    // Loop over random points within the cell
    for (int irnd = 0; irnd < AtlasGeoTestsConfig::N_RANDOM_POINTS; ++irnd) {
      Hit2D hit {};
      hit.m_eta = AtlasGeoTestsConfig::sample(cell.eta() - 0.49 * cell.deta(),
                                              cell.eta() + 0.49 * cell.deta());
      hit.m_phi = Cell::norm_angle(AtlasGeoTestsConfig::sample(
          cell.phi() - 0.49 * cell.dphi(), cell.phi() + 0.49 * cell.dphi()));

      // Retrieve the cell closest to the hit position
      const auto& best_cell = AtlasGeoTests::geo->get_cell(layer, hit);

      // Check if the hit position is inside the cell
      ASSERT_NEAR(hit.m_eta, best_cell.eta(), 0.5 * best_cell.deta());
      ASSERT_NEAR(Cell::norm_angle(hit.m_phi - best_cell.phi()),
                  0.0,
                  0.5 * best_cell.dphi());
    }
  }
}

TEST_F(AtlasGeoTests, FCAL0Lookup)

{
  constexpr int layer = 21;

  // Loop over all cells in the layer
  for (int i = 0; i < AtlasGeoTests::geo->n_cells(layer); ++i) {
    const auto& cell = AtlasGeoTests::geo->get_cell_at_idx(layer, i);

    // Loop over random points within the cell
    for (int irnd = 0; irnd < AtlasGeoTestsConfig::N_RANDOM_POINTS; ++irnd) {
      Hit2D hit {};
      hit.m_x = AtlasGeoTestsConfig::sample(cell.x() - 0.49 * cell.dx(),
                                            cell.x() + 0.49 * cell.dx());
      hit.m_y = AtlasGeoTestsConfig::sample(cell.y() - 0.49 * cell.dy(),
                                            cell.y() + 0.49 * cell.dy());

      // Retrieve the cell closest to the hit position
      const auto& best_cell = AtlasGeoTests::geo->get_cell(layer, hit);

      // Check if the hit position is inside the cell
      ASSERT_NEAR(hit.m_x, best_cell.x(), 0.5 * best_cell.dx());
      ASSERT_NEAR(hit.m_y, best_cell.y(), 0.5 * best_cell.dy());
    }
  }
}
