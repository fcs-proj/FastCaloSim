
#include <CLHEP/Random/RanluxEngine.h>
#include <gtest/gtest.h>

#include "AtlasGeoTests.h"
#include "TFile.h"
#include "TestConfig/BasicSimTestsConfig.h"

class BasicSimTests : public AtlasGeoTests
{
protected:
  static TFile* param_file;

  // Per-test-suite set-up.
  // Called before the first test in this test suite.
  // Note: this does not work if incurred via ctest
  static void SetUpTestSuite()
  {
    // Call the base setup that sets up the geometry
    AtlasGeoTests::SetUpTestSuite();

    // Load parametrization file
    param_file =
        new TFile((BasicSimTestsConfig::PARAM_FILE_PATH).c_str(), "READ");
  }

  // Per-test-suite tear-down.
  // Called after the last test in this test suite.
  // Note: this does not work if incurred via ctest
  static void TearDownTestSuite()
  {
    // Call the base tear down that deletes the geometry
    AtlasGeoTests::TearDownTestSuite();

    param_file->Close();
    delete param_file;
  }
};

// Initialize the static members
TFile* BasicSimTests::param_file = nullptr;