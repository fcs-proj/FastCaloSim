// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#include <memory>
#include <string>

#include <TClass.h>
#include <TMemFile.h>
#include <gtest/gtest.h>

#include "FastCaloSim/Core/TFCSBinnedShower.h"
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

// With an HDF5 file path, the event library is not written to the file but
// stays in memory, and reading the object back must not need the geometry
// (it is set only afterwards; reading used to dereference it)
TEST(BinnedShowerTests, EventLibraryFromHDF5IsNotWritten)
{
  const std::unique_ptr<TFCSParametrizationBase> param(
      static_cast<TFCSParametrizationBase*>(
          TClass::GetClass("TFCSBinnedShower")->New()));
  ASSERT_NE(param, nullptr);
  auto* shower = static_cast<TFCSBinnedShower*>(param.get());
  shower->set_event_library(TFCSBinnedShower::eventvector_t(3));
  const std::string path = "event_library.h5";
  shower->set_hdf5_path(path);

  TMemFile file("binned_shower.root", "RECREATE");
  ASSERT_GT(file.WriteTObject(param.get(), "shower"), 0);
  EXPECT_EQ(shower->get_eventlibrary().size(), 3U);

  const std::unique_ptr<TFCSParametrizationBase> read(
      file.Get<TFCSParametrizationBase>("shower"));
  ASSERT_NE(read, nullptr);
  const auto* read_shower = static_cast<const TFCSBinnedShower*>(read.get());
  EXPECT_TRUE(read_shower->get_eventlibrary().empty());
  EXPECT_EQ(read_shower->get_hdf5_path(), path);
}
