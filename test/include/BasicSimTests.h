// Copyright (c) 2024 CERN for the benefit of the FastCaloSim project

#include <CLHEP/Random/RanluxEngine.h>
#include <gtest/gtest.h>

#include "AtlasGeoTests.h"
#include "TFile.h"
#include "TestConfig/BasicSimTestsConfig.h"

class BasicSimTests : public AtlasGeoTests
{
protected:
  static std::map<std::string, TFile*> param_files;

  // Per-test-suite set-up.
  // Called before the first test in this test suite.
  static void SetUpTestSuite()
  {
    // Call the base setup that sets up the geometry
    AtlasGeoTests::SetUpTestSuite();

    // Load the parametrization files for the different regions
    for (const auto& [region, path] : BasicSimTestsConfig::PARAM_FILE_PATHS) {
      param_files[region] = new TFile(path.c_str(), "READ");
    }
  }

  // Per-test-suite tear-down.
  // Called after the last test in this test suite.
  static void TearDownTestSuite()
  {
    // Call the base tear down that deletes the geometry
    AtlasGeoTests::TearDownTestSuite();

    // Close and delete the parametrization files
    for (const auto& [region, file] : param_files) {
      file->Close();
      delete file;
    }
  }
};

// Initialize the static members
std::map<std::string, TFile*> BasicSimTests::param_files;
