// Copyright (c) 2024 CERN for the benefit of the FastCaloSim project

#pragma once

#include "BasicExtrapolTests.h"
#include "BasicSimTests.h"
#include "G4RunTests.h"

class AtlasSimTests
    : public BasicExtrapolTests
    , public BasicSimTests
    , public G4RunTests
{
protected:
  static void SetUpTestSuite()
  {
    // Call setup for G4RunTests which gives us the basic Geant4 environment
    G4RunTests::SetUpTestSuite();

    // Call setup for BasicSimTests which gives us the set up parametrization
    BasicSimTests::SetUpTestSuite();

    // Call setup for BasicExtrapolTests which gives us the set up extrapolation
    BasicExtrapolTests::SetUpTestSuite();
  }

  static void TearDownTestSuite()
  {
    // Call teardown for G4RunTests
    G4RunTests::TearDownTestSuite();

    // Call teardown for BasicSimTests
    BasicSimTests::TearDownTestSuite();

    // Call teardown for BasicExtrapolTests
    BasicExtrapolTests::TearDownTestSuite();
  }
};
