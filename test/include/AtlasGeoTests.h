// Copyright (c) 2024 CERN for the benefit of the FastCaloSim project

#pragma once

#include <filesystem>

#include <CLHEP/Random/RanluxEngine.h>
#include <gtest/gtest.h>
#include <unistd.h>

#include "ATLASFCalGeoPlugin/FCal.h"
#include "FastCaloSim/Geometry/CaloGeo.h"
#include "ROOT/RDataFrame.hxx"
#include "TestConfig/AtlasGeoTestsConfig.h"
#include "TestHelpers/Event.h"

// Test suite for GeoTests
class AtlasGeoTests
    : public virtual ::testing::TestWithParam<TestHelpers::Event>
{
protected:
  static CaloGeo* geo;

  // Sets up the test suite
  // Called before the first test in this test suite
  static void SetUpTestSuite()
  {
    ROOT::RDataFrame df =
        ROOT::RDataFrame(AtlasGeoTestsConfig::ATLAS_CALO_CELL_TREE_NAME,
                         AtlasGeoTestsConfig::ATLAS_CALO_CELL_PATH);
    geo = new CaloGeo();

    // Create Temporary dir for the RTree files
    char tmpl[] = "/tmp/fastcalosimXXXXXX";
    int fd = mkstemp(tmpl);
    close(fd);
    std::remove(tmpl);
    std::filesystem::create_directory(tmpl);
    std::filesystem::path tmp_dir_path {tmpl};

    // Build the geometry and write it to disk
    geo->build(df, tmp_dir_path.string());
    // Load the geometry from the files
    // Set the Rtree cache size to 5 MB per layer
    size_t rtree_cache_size = 5 * 1024 * 1024;
    // Set the cell store cache size to 10 MB
    size_t cell_cache_size = 10 * 1024 * 1024;
    geo->load(tmp_dir_path.string(), rtree_cache_size, cell_cache_size);

    // Create alternative geometry handler for FCal
    std::shared_ptr<FCal> fcal_geo = std::make_shared<FCal>();
    // Load the FCal geometry from the files
    fcal_geo->load(AtlasGeoTestsConfig::FCAL_ELECTRODE_FILES);
    // Set the pointer to the main geo to access cell store
    fcal_geo->set_geo(geo);
    // Set the alternative geometry handler for the FCal layers (21 - 23)
    geo->set_alt_geo_handler(21, 23, fcal_geo);
  }

  // Tears down the test suite
  // Called after the last test in this test suite
  static void TearDownTestSuite() {}  // deleting geo here causes issues
};

// Initialize the static member
CaloGeo* AtlasGeoTests::geo = nullptr;
