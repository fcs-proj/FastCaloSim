// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#include <stdexcept>

#include "FastCaloSim/Definitions/ParticleData.h"

#include <gtest/gtest.h>

// Charges of the built-in default provider, used when no experiment provider
// is installed
TEST(ParticleDataTests, DefaultProviderCharges)
{
  ParticleData::setProvider(nullptr);

  EXPECT_EQ(ParticleData::charge(11), -1.);  // electron
  EXPECT_EQ(ParticleData::charge(-11), 1.);  // positron
  EXPECT_EQ(ParticleData::charge(211), 1.);  // pi+
  EXPECT_EQ(ParticleData::charge(-211), -1.);
  EXPECT_EQ(ParticleData::charge(321), 1.);  // K+
  EXPECT_EQ(ParticleData::charge(-321), -1.);
  EXPECT_EQ(ParticleData::charge(2212), 1.);  // proton
  EXPECT_EQ(ParticleData::charge(-2212), -1.);

  for (const int neutral : {22, 2112, -2112, 111, 130, 310}) {
    EXPECT_EQ(ParticleData::charge(neutral), 0.) << "pdgID " << neutral;
  }

  EXPECT_THROW((void)ParticleData::charge(13), std::runtime_error);
}
