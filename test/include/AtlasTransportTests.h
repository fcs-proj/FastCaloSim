#pragma once

#include "BasicExtrapolTests.h"
#include "BasicSimTests.h"
#include "G4RunTests.h"

class AtlasTransportTests
    : public G4RunTests
    , public BasicExtrapolTests
{
protected:
  static void SetUpTestSuite()
  {
    // Call setup for G4RunTests which gives us the basic Geant4 environment
    G4RunTests::SetUpTestSuite();

    // Call setup for BasicExtrapolTests which gives us the set up extrapolation
    BasicExtrapolTests::SetUpTestSuite();
  }

  static void TearDownTestSuite()
  {
    // Call teardown for G4RunTests
    G4RunTests::TearDownTestSuite();

    // Call teardown for BasicExtrapolTests
    BasicExtrapolTests::TearDownTestSuite();
  }
};