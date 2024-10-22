// Copyright (c) 2024 CERN for the benefit of the FastCaloSim project

#pragma once

#include <CLHEP/Random/RanluxEngine.h>
#include <gtest/gtest.h>

#include "TestHelpers/Event.h"

class TFCSSimulationStateTest
    : public ::testing::TestWithParam<TestHelpers::Event>
{
protected:
  CLHEP::HepRandomEngine* m_random_engine;

  void SetUp() override { m_random_engine = new CLHEP::RanluxEngine(); }

  void TearDown() override { delete m_random_engine; }
};
