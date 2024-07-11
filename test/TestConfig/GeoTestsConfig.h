#pragma once

#include <gtest/gtest.h>

class GeoTestsConfig : public ::testing::Environment
{
public:
  // The location of the ATLAS calo geometry file
  inline static const std::string ATLAS_CALO_CELL_PATH =
      std::string(TEST_BASE_DIR) + "/data/geometry/ATLASCaloCells.root";
  // The name of the tree in the ATLAS calo geometry file
  inline static const std::string ATLAS_CALO_CELL_TREE_NAME = "caloDetCells";
};
