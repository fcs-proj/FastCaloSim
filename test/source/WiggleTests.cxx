// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#include <memory>

#include <CLHEP/Random/MixMaxRng.h>
#include <TClass.h>
#include <gtest/gtest.h>

#include "FastCaloSim/Core/TFCSHitCellMappingWiggle.h"
#include "FastCaloSim/Core/TFCSParametrizationBase.h"
#include "FastCaloSim/Core/TFCSSimulationState.h"
#include "FastCaloSim/Geometry/CaloGeo.h"
#include "FastCaloSim/Geometry/Cell.h"

namespace
{
// A geometry with one cell around eta = phi = 0
class OneCellGeo : public CaloGeo
{
public:
  auto get_cell(std::size_t /*layer*/, const Position& /*pos*/) const
      -> Cell override
  {
    const Position center {
        .m_x = 0, .m_y = 0, .m_z = 0, .m_eta = 0, .m_phi = 0, .m_r = 0};
    return {7,
            center,
            0,
            /*isBarrel=*/true,
            /*isXYZ=*/false,
            /*isEtaPhiR=*/true,
            /*isEtaPhiZ=*/false,
            /*isRPhiZ=*/false,
            0,
            0,
            0,
            1.0,
            1.0,
            0};
  }
};
}  // namespace

// After clear() there are no wiggle bins: a hit is mapped to its cell without
// a wiggle, instead of reading bins that do not exist
TEST(WiggleTest, MapsHitsAfterClear)
{
  // Created through ROOT, as when read from a parametrization file. The class
  // is not exported, so the test owns it through the exported base class.
  const std::unique_ptr<TFCSParametrizationBase> owner(
      static_cast<TFCSParametrizationBase*>(
          TClass::GetClass("TFCSHitCellMappingWiggle")->New()));
  auto& wiggle = *static_cast<TFCSHitCellMappingWiggle*>(owner.get());
  // A new mapping has one bin without a function
  ASSERT_EQ(wiggle.get_number_of_bins(), 1U);
  wiggle.clear();
  ASSERT_EQ(wiggle.get_number_of_bins(), 0U);
  OneCellGeo geo;
  owner->set_geometry(&geo);

  CLHEP::MixMaxRng engine(42);
  TFCSSimulationState simulstate(&engine);
  TFCSLateralShapeParametrizationHitBase::Hit hit(0.1, 0.2, 5.0);

  // Called through the base class, as the hit chains do
  TFCSLateralShapeParametrizationHitBase& mapping = wiggle;
  EXPECT_EQ(mapping.simulate_hit(hit, simulstate, nullptr, nullptr),
            FCSSuccess);
  EXPECT_DOUBLE_EQ(hit.phi(), 0.2);
  ASSERT_EQ(simulstate.cells().size(), 1U);
  EXPECT_DOUBLE_EQ(simulstate.cells().at(7), 5.0);
}
