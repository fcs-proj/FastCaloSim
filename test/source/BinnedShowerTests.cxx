// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#include <memory>

#include <TClass.h>
#include <gtest/gtest.h>

#include "FastCaloSim/Core/TFCSExtrapolationState.h"
#include "FastCaloSim/Core/TFCSParametrizationBase.h"
#include "FastCaloSim/Core/TFCSSimulationState.h"
#include "FastCaloSim/Core/TFCSTruthState.h"
#include "FastCaloSim/Geometry/CaloGeo.h"

// TFCSBinnedShower::simulate() first removes the data of the previous event
// from the simulation state. On a new state there is none, which used to throw
// std::out_of_range.
TEST(BinnedShowerTests, SimulateOnNewState)
{
  // Created through the dictionary and owned through the exported base class
  const std::unique_ptr<TFCSParametrizationBase> param(
      static_cast<TFCSParametrizationBase*>(
          TClass::GetClass("TFCSBinnedShower")->New()));
  ASSERT_NE(param, nullptr);
  CaloGeo geo;  // no layers
  param->set_geometry(&geo);

  TFCSSimulationState simulstate;
  const TFCSTruthState truth;
  const TFCSExtrapolationState extrapol;
  // The event library is empty: the simulation reports an error and deposits
  // nothing
  EXPECT_NO_THROW(param->simulate(simulstate, &truth, &extrapol));
  EXPECT_DOUBLE_EQ(simulstate.E(), 0.0);
}
