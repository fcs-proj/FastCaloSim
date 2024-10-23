// Copyright (c) 2024 CERN for the benefit of the FastCaloSim project

#pragma once

#include <TestHelpers/Event.h>
#include <gtest/gtest.h>

class TFCSSTruthStateTest : public ::testing::TestWithParam<TestHelpers::Event>
{
protected:
  void SetUp() override {}

  void TearDown() override {}
};
