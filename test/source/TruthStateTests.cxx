// Copyright (c) 2024 CERN for the benefit of the FastCaloSim project

#include "TruthStateTests.h"

#include <gtest/gtest.h>

#include "FastCaloSim/Core/TFCSTruthState.h"

// Test default constructor
TEST(TFCSTruthStateTest, DefaultConstructor)
{
  TFCSTruthState ts;
  EXPECT_EQ(ts.pdgid(), 0);
  EXPECT_EQ(ts.vertex(), TLorentzVector(0, 0, 0, 0));
  EXPECT_DOUBLE_EQ(ts.Ekin_off(), 0.0);
}

// Test parameterized constructor
TEST(TFCSTruthStateTest, ParamConstructor)
{
  TFCSTruthState ts(1.0, 2.0, 3.0, 4.0, 5);
  EXPECT_EQ(ts.pdgid(), 5);
  EXPECT_EQ(ts.X(), 1.0);
  EXPECT_EQ(ts.Y(), 2.0);
  EXPECT_EQ(ts.Z(), 3.0);
  EXPECT_EQ(ts.T(), 4.0);
}

// Test setters and getters
TEST(TFCSTruthStateTest, SettersAndGetters)
{
  TFCSTruthState ts;
  ts.set_pdgid(10);
  EXPECT_EQ(ts.pdgid(), 10);

  ts.set_vertex(1.0, 2.0, 3.0, 4.0);
  EXPECT_EQ(ts.vertex(), TLorentzVector(1.0, 2.0, 3.0, 4.0));

  ts.set_Ekin_off(5.0);
  EXPECT_DOUBLE_EQ(ts.Ekin_off(), 5.0);
}

// Test Ekin computation
TEST(TFCSTruthStateTest, EkinComputation)
{
  TFCSTruthState ts(0, 0, 1, sqrt(1 + 0.25), 0);  // mass = 0.5
  ts.set_Ekin_off(0.1);
  double expectedEkin = sqrt(1 + 0.25) - 0.5 + 0.1;  // E - M + m_ekin_off
  EXPECT_NEAR(ts.Ekin(), expectedEkin, 1e-6);
}

// Print functionality (mocking FCS_MSG_INFO)
TEST(TFCSTruthStateTest, Print)
{
  TFCSTruthState ts(1.0, 0, 0, 1.0, 11);
  testing::internal::CaptureStdout();
  ts.Print("");
  std::string output = testing::internal::GetCapturedStdout();
  EXPECT_TRUE(output.find("PDGID=11") != std::string::npos);
}
