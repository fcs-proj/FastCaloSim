// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#include <gtest/gtest.h>

#include "FastCaloSim/Core/TFCSExtrapolationState.h"
#include "FastCaloSim/Geometry/SubPos.h"

using FastCaloSim::SubPos;

// Positions the extrapolation did not reach read as the defaults
TEST(ExtrapolationStateTests, DefaultsForPositionsNotSet)
{
  const TFCSExtrapolationState state;
  EXPECT_FALSE(state.OK(2, SubPos::MID));
  EXPECT_DOUBLE_EQ(state.eta(2, SubPos::MID), -999);
  EXPECT_DOUBLE_EQ(state.phi(2, SubPos::ENT), -999);
  EXPECT_DOUBLE_EQ(state.r(2, SubPos::EXT), 0);
  EXPECT_DOUBLE_EQ(state.z(2, SubPos::MID), 0);
  EXPECT_DOUBLE_EQ(state.d(2, SubPos::MID), 0);
  EXPECT_DOUBLE_EQ(state.detaBorder(2, SubPos::MID), 0);
}

// Each value of a position is set on its own; layers below the highest set
// layer and the other sub-positions keep the defaults
TEST(ExtrapolationStateTests, SetValues)
{
  TFCSExtrapolationState state;
  state.set_eta(5, SubPos::ENT, 1.5);
  state.set_OK(5, SubPos::ENT);
  state.set_r(5, SubPos::EXT, 1500);

  EXPECT_DOUBLE_EQ(state.eta(5, SubPos::ENT), 1.5);
  EXPECT_TRUE(state.OK(5, SubPos::ENT));
  EXPECT_DOUBLE_EQ(state.phi(5, SubPos::ENT), -999);
  EXPECT_DOUBLE_EQ(state.r(5, SubPos::EXT), 1500);
  EXPECT_FALSE(state.OK(5, SubPos::EXT));
  EXPECT_DOUBLE_EQ(state.eta(5, SubPos::MID), -999);
  EXPECT_DOUBLE_EQ(state.eta(3, SubPos::ENT), -999);
  EXPECT_DOUBLE_EQ(state.eta(23, SubPos::ENT), -999);
}

TEST(ExtrapolationStateTests, ClearResetsAll)
{
  TFCSExtrapolationState state;
  state.set_eta(1, SubPos::MID, 0.5);
  state.set_OK(1, SubPos::MID);
  state.set_IDCaloBoundary_eta(0.4);
  state.clear();
  EXPECT_FALSE(state.OK(1, SubPos::MID));
  EXPECT_DOUBLE_EQ(state.eta(1, SubPos::MID), -999);
  EXPECT_DOUBLE_EQ(state.IDCaloBoundary_eta(), -999);
}
