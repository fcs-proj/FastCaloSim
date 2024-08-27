#pragma once

#include <TestHelpers/Event.h>
#include <gtest/gtest.h>

class TFCSSTruthStateTest : public ::testing::TestWithParam<TestHelpers::Event>
{
protected:
  void SetUp() override {}

  void TearDown() override {}
};