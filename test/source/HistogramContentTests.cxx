// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#include <cstdint>
#include <limits>

#include <gtest/gtest.h>

#include "FastCaloSim/Core/TFCS1DFunctionTemplateHelpers.h"

// A fraction of 1 is stored as the largest value of the type, also for 32 bit,
// where the largest value rounds up to 2^32 in single precision
TEST(HistogramContentTest, FractionOneSaturates)
{
  using Numeric32 = TFCS1DFunction_Numeric<std::uint32_t, float>;
  EXPECT_EQ(Numeric32::MaxCeilOnlyForInt(1.0F),
            std::numeric_limits<std::uint32_t>::max());
  EXPECT_EQ(Numeric32::MaxCeilOnlyForInt(0.0F), 0U);
  EXPECT_EQ(Numeric32::MaxCeilOnlyForInt(0.5F), 2147483648U);

  // Unchanged for the smaller types
  EXPECT_EQ(
      (TFCS1DFunction_Numeric<std::uint8_t, float>::MaxCeilOnlyForInt(1.0F)),
      255U);
  EXPECT_EQ(
      (TFCS1DFunction_Numeric<std::uint16_t, float>::MaxCeilOnlyForInt(1.0F)),
      65535U);

  TFCS1DFunction_HistogramContent<std::uint32_t, float> content;
  content.set_nbins(3);
  content.set_fraction(0, 0.25F);
  content.set_fraction(1, 1.0F);
  EXPECT_FLOAT_EQ(content.get_fraction(0), 0.25F);
  EXPECT_FLOAT_EQ(content.get_fraction(1), 1.0F);
}
