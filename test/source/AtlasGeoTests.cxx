#include "AtlasGeoTests.h"

#include <gtest/gtest.h>

#include "FastCaloSim/Geometry/CaloGeo.h"

struct Hit
{
  auto x() const -> double { return m_x; }
  auto y() const -> double { return m_y; }
  auto z() const -> double { return m_z; }
  auto eta() const -> double { return m_eta; }
  auto phi() const -> double { return m_phi; }
  auto r() const -> double { return m_r; }
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
      Hit hit {};
      hit.m_eta = AtlasGeoTestsConfig::sample(cell.eta() - 0.49 * cell.deta(),
                                              cell.eta() + 0.49 * cell.deta());
      hit.m_phi = AtlasGeoTestsConfig::normalize_angle(
          AtlasGeoTestsConfig::sample(cell.phi() - 0.49 * cell.dphi(),
                                      cell.phi() + 0.49 * cell.dphi()));
      hit.m_r = AtlasGeoTestsConfig::sample(cell.r() - 0.49 * cell.dr(),
                                            cell.r() + 0.49 * cell.dr());

      // Retrieve the cell closest to the hit position
      const auto& best_cell = AtlasGeoTests::geo->get_cell(layer, hit);

      // Check if the hit position is inside the cell
      ASSERT_NEAR(hit.m_eta, best_cell.eta(), 0.5 * best_cell.deta());
      ASSERT_NEAR(
          AtlasGeoTestsConfig::normalize_angle(hit.m_phi - best_cell.phi()),
          0.0,
          0.5 * best_cell.dphi());
      ASSERT_NEAR(hit.m_r, best_cell.r(), 0.5 * best_cell.dr());
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
      Hit hit {};
      hit.m_eta = AtlasGeoTestsConfig::sample(cell.eta() - 0.49 * cell.deta(),
                                              cell.eta() + 0.49 * cell.deta());
      hit.m_phi = AtlasGeoTestsConfig::normalize_angle(
          AtlasGeoTestsConfig::sample(cell.phi() - 0.49 * cell.dphi(),
                                      cell.phi() + 0.49 * cell.dphi()));
      hit.m_z = AtlasGeoTestsConfig::sample(cell.z() - 0.49 * cell.dz(),
                                            cell.z() + 0.49 * cell.dz());

      // Retrieve the cell closest to the hit position
      const auto& best_cell = AtlasGeoTests::geo->get_cell(layer, hit);

      // Check if the hit position is inside the cell
      ASSERT_NEAR(hit.m_eta, best_cell.eta(), 0.5 * best_cell.deta());
      ASSERT_NEAR(
          AtlasGeoTestsConfig::normalize_angle(hit.m_phi - best_cell.phi()),
          0.0,
          0.5 * best_cell.dphi());
      ASSERT_NEAR(hit.m_z, best_cell.z(), 0.5 * best_cell.dz());
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
      Hit hit {};
      hit.m_x = AtlasGeoTestsConfig::sample(cell.x() - 0.49 * cell.dx(),
                                            cell.x() + 0.49 * cell.dx());
      hit.m_y = AtlasGeoTestsConfig::sample(cell.y() - 0.49 * cell.dy(),
                                            cell.y() + 0.49 * cell.dy());
      hit.m_z = AtlasGeoTestsConfig::sample(cell.z() - 0.49 * cell.dz(),
                                            cell.z() + 0.49 * cell.dz());

      // Retrieve the cell closest to the hit position
      const auto& best_cell = AtlasGeoTests::geo->get_cell(layer, hit);

      // Check if the hit position is inside the cell
      ASSERT_NEAR(hit.m_x, best_cell.x(), 0.5 * best_cell.dx());
      ASSERT_NEAR(hit.m_y, best_cell.y(), 0.5 * best_cell.dy());
      ASSERT_NEAR(hit.m_z, best_cell.z(), 0.5 * best_cell.dz());
    }
  }
}
