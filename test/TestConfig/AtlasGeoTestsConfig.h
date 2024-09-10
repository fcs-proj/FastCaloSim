#pragma once

#include <TRandom.h>
#include <gtest/gtest.h>

class AtlasGeoTestsConfig : public ::testing::Environment
{
public:
  /// @brief The location of the ATLAS calo geometry file
  inline static const std::string ATLAS_CALO_CELL_PATH =
      std::string(TEST_BASE_DIR) + "/data/geometry/ATLAS/CaloCells.root";

  /// @brief The location of the ATLAS FCal files
  /// Note: as the FCal cells are irregular and cannot be described by simple
  /// cuboids for the hit-cell matching, the geometry is stored in separate
  /// files
  inline static const std::string FCAL1_GEO_PATH = std::string(TEST_BASE_DIR)
      + "/data/geometry/ATLAS/FCal/FCal1-electrodes.sorted.HV.09Nov2007.dat";
  inline static const std::string FCAL2_GEO_PATH = std::string(TEST_BASE_DIR)
      + "/data/geometry/ATLAS/FCal/FCal2-electrodes.sorted.HV.April2011.dat";
  inline static const std::string FCAL3_GEO_PATH = std::string(TEST_BASE_DIR)
      + "/data/geometry/ATLAS/FCal/FCal3-electrodes.sorted.HV.09Nov2007.dat";

  inline static const std::array<std::string, 3> FCAL_ELECTRODE_FILES = {
      FCAL1_GEO_PATH, FCAL2_GEO_PATH, FCAL3_GEO_PATH};

  // The name of the tree in the ATLAS calo geometry file
  inline static const std::string ATLAS_CALO_CELL_TREE_NAME = "caloDetCells";
  // The number of random points that will be generated within a cell for
  // testing lookup
  inline static const int N_RANDOM_POINTS = 100;

  // Helper method to sample a random number between a and b
  static auto sample(double a, double b) -> double
  {
    return a + (b - a) * gRandom->Rndm();
  }
};
