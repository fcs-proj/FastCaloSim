#include <string>

#include <gtest/gtest.h>

#include "FastCaloSim/FastCaloSim.h"

TEST(BasicClassTests, ClassExportOk)
{
  auto const exported = exported_class {};

  auto const expected_name = std::string("FastCaloSim");
  const std::string actual_name = exported.name();

  EXPECT_EQ(expected_name, actual_name);
}
