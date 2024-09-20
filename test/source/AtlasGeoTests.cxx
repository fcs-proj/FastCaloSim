#include "AtlasGeoTests.h"

#include <gtest/gtest.h>

#include "FastCaloSim/Geometry/CaloGeo.h"

TEST_F(AtlasGeoTests, PreSamplerBLookup)
{
  constexpr int layer = 0;

  // Loop over all cells in the layer
  for (int i = 0; i < AtlasGeoTests::geo->n_cells(layer); ++i) {
    const auto& cell = AtlasGeoTests::geo->get_cell_at_idx(layer, i);

    // Loop over random points within the cell
    for (int irnd = 0; irnd < AtlasGeoTestsConfig::N_RANDOM_POINTS; ++irnd) {
      Position pos {};
      pos.m_eta = AtlasGeoTestsConfig::sample(cell.eta() - 0.49 * cell.deta(),
                                              cell.eta() + 0.49 * cell.deta());
      pos.m_phi = Cell::norm_angle(AtlasGeoTestsConfig::sample(
          cell.phi() - 0.49 * cell.dphi(), cell.phi() + 0.49 * cell.dphi()));

      // Retrieve the cell closest to the hit position
      const auto& best_cell = AtlasGeoTests::geo->get_cell(layer, pos);

      // Check if the hit position is inside the cell
      ASSERT_NEAR(pos.m_eta, best_cell.eta(), 0.5 * best_cell.deta());
      ASSERT_NEAR(Cell::norm_angle(pos.m_phi - best_cell.phi()),
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
      Position pos {};
      pos.m_eta = AtlasGeoTestsConfig::sample(cell.eta() - 0.49 * cell.deta(),
                                              cell.eta() + 0.49 * cell.deta());
      pos.m_phi = Cell::norm_angle(AtlasGeoTestsConfig::sample(
          cell.phi() - 0.49 * cell.dphi(), cell.phi() + 0.49 * cell.dphi()));

      // Retrieve the cell closest to the hit position
      const auto& best_cell = AtlasGeoTests::geo->get_cell(layer, pos);

      // Check if the hit position is inside the cell
      ASSERT_NEAR(pos.m_eta, best_cell.eta(), 0.5 * best_cell.deta());
      ASSERT_NEAR(Cell::norm_angle(pos.m_phi - best_cell.phi()),
                  0.0,
                  0.5 * best_cell.dphi());
    }
  }
}
