// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#include <cstddef>
#include <stdexcept>

#include "AtlasGeoTests.h"

#include <gtest/gtest.h>

#include "FastCaloSim/Geometry/CaloGeo.h"

TEST_F(AtlasGeoTests, PreSamplerBLookup)
{
  constexpr int layer = 0;

  // Loop over all cells in the layer
  for (unsigned int i = 0; i < AtlasGeoTests::geo->n_cells(layer); ++i) {
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
  for (unsigned int i = 0; i < AtlasGeoTests::geo->n_cells(layer); ++i) {
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

TEST_F(AtlasGeoTests, LayerCoordinateSystems)
{
  // The parametrizations simulate showers in (x, y) in the layers whose cells
  // are described in (x, y, z), which are the FCal layers 21 to 23 in ATLAS
  constexpr std::size_t n_layers = 24;
  ASSERT_EQ(AtlasGeoTests::geo->n_layers(), n_layers);
  for (std::size_t layer = 0; layer < n_layers; ++layer) {
    const bool xyz = AtlasGeoTests::geo->is_xyz(layer);
    EXPECT_EQ(xyz, layer >= 21) << "layer " << layer;
    EXPECT_EQ(static_cast<int>(xyz)
                  + static_cast<int>(AtlasGeoTests::geo->is_eta_phi_r(layer))
                  + static_cast<int>(AtlasGeoTests::geo->is_eta_phi_z(layer)),
              1)
        << "layer " << layer;
  }
  EXPECT_THROW(static_cast<void>(AtlasGeoTests::geo->is_xyz(n_layers)),
               std::out_of_range);
}

TEST_F(AtlasGeoTests, FCalPositionWithoutCell)
{
  constexpr std::size_t layer = 21;

  // A position in an FCal cell is found by both lookups
  const Cell cell = AtlasGeoTests::geo->get_cell_at_idx(layer, 0);
  const Position in_cell {.m_x = cell.x(),
                          .m_y = cell.y(),
                          .m_z = cell.z(),
                          .m_eta = 0,
                          .m_phi = 0,
                          .m_r = 0};
  EXPECT_EQ(AtlasGeoTests::geo->find_cell(layer, in_cell).id(), cell.id());
  EXPECT_EQ(AtlasGeoTests::geo->get_cell(layer, in_cell).id(), cell.id());

  // In the beam hole there is no FCal cell: find_cell gives an invalid cell,
  // so that TFCSHitCellMappingFCal can simulate the hit again, and get_cell
  // throws
  const Position in_beam_hole {
      .m_x = 1, .m_y = 1, .m_z = cell.z(), .m_eta = 0, .m_phi = 0, .m_r = 0};
  EXPECT_FALSE(AtlasGeoTests::geo->find_cell(layer, in_beam_hole).is_valid());
  EXPECT_THROW(AtlasGeoTests::geo->get_cell(layer, in_beam_hole),
               std::runtime_error);
}
