// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#include <TInterpreter.h>
#include <gtest/gtest.h>

// The ROOT interpreter parses the FastCaloSim headers when an interactive
// session or macro uses a class: the headers must be accepted next to the
// forward declarations that the dictionary provides.
TEST(RootInterpreterTest, TemplateHistogramHelpers)
{
  EXPECT_EQ(gInterpreter->ProcessLine(
                "TFCS1DFunction_Numeric<unsigned char, float>::MaxValue;"),
            255);
}
