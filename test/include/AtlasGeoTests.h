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

    // 5mb cache size per layer for the RTree
    geo->build(df, tmp_dir_path.string(), true, 5 * 1024 * 1024);

    // Create alternative geometry handler for FCal
    std::shared_ptr<FCal> fcal_geo = std::make_shared<FCal>();
    // Load the FCal geometry from the files
    fcal_geo->load(AtlasGeoTestsConfig::FCAL_ELECTRODE_FILES);
    // Cast the FCal geometry handler to the geo interface
    std::shared_ptr<CaloGeo> fcal_geo_handler =
        std::static_pointer_cast<CaloGeo>(fcal_geo);
    // Set the alternative geometry handler for the FCal layers (21 - 23)
    geo->set_alt_geo_handler(21, 23, fcal_geo_handler);
  }

  // Tears down the test suite
  // Called after the last test in this test suite
  static void TearDownTestSuite() {}  // deleting geo here causes issues
};

// Initialize the static member
CaloGeo* AtlasGeoTests::geo = nullptr;
