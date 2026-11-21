// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#include <TInterpreter.h>
#include <TSystem.h>
#include <gtest/gtest.h>

// The ROOT interpreter parses the FastCaloSim headers when an interactive
// session or macro uses a class: the headers must be accepted next to the
// forward declarations that the dictionary provides. The library is loaded
// explicitly, as in an interactive session: the test executable uses none of
// its symbols, so the linker may drop it (--as-needed).
TEST(RootInterpreterTest, TemplateHistogramHelpers)
{
  ASSERT_GE(gSystem->Load(FastCaloSimLib), 0);
  EXPECT_EQ(gInterpreter->ProcessLine(
                "TFCS1DFunction_Numeric<unsigned char, float>::MaxValue;"),
            255);
}
