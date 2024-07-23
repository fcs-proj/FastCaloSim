#include "AtlasDeprecatedGeoTests.h"

#include <CLHEP/Random/RanluxEngine.h>
#include <gtest/gtest.h>

TEST_F(AtlasDeprecatedGeoTests, GetGeoCellBySamplingEtaPhi)
{
  const CaloDetDescrElement* cell = geo->getDDE(2, 0.24, 0.24);
  long int cell_id = cell->identify();
  EXPECT_EQ(cell_id, 3296654795753914368);
  EXPECT_EQ(cell->getSampling(), 2);
  EXPECT_NEAR(cell->eta(), 0.238859, 1e-4);
  EXPECT_NEAR(cell->phi(), 0.228048, 1e-4);
}

TEST_F(AtlasDeprecatedGeoTests, GetGeoCellByIdentifier)
{
  const CaloDetDescrElement* cell = geo->getDDE(3260641881524011008);
  long int cell_id = cell->identify();
  EXPECT_EQ(cell_id, 3260641881524011008);
  EXPECT_EQ(cell->getSampling(), 2);
  EXPECT_NEAR(cell->eta(), -0.411563, 1e-4);
  EXPECT_NEAR(cell->phi(), 1.63532, 1e-4);
}
