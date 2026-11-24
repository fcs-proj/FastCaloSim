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
#include "FastCaloSim/Geometry/LayerRoles.h"

namespace
{
// Simulates with a binned shower of the given class that has no shower to take
// an event from (no event library, or no model, loaded)
auto simulate_without_showers(const char* class_name, double& energy)
    -> FCSReturnCode
{
  // Created through the dictionary and owned through the exported base class
  const std::unique_ptr<TFCSParametrizationBase> param(
      static_cast<TFCSParametrizationBase*>(
          TClass::GetClass(class_name)->New()));
  EXPECT_NE(param, nullptr);
  CaloGeo geo;  // no layers
  FastCaloSim::LayerRoles roles;
  roles.shower_center_layers = {{.layer = 0, .max_abs_eta = 0}};
  geo.set_layer_roles(roles);
  param->set_geometry(&geo);

  TFCSSimulationState simulstate;
  const TFCSTruthState truth;
  const TFCSExtrapolationState extrapol;
  const FCSReturnCode status = param->simulate(simulstate, &truth, &extrapol);
  energy = simulstate.E();
  return status;
}
}  // namespace

// TFCSBinnedShower::simulate() first removes the data of the previous event
// from the simulation state. On a new state there is none, which used to throw
// std::out_of_range. The event library is empty: the simulation fails and
// deposits nothing.
TEST(BinnedShowerTests, SimulateOnNewState)
{
  double energy = -1;
  EXPECT_EQ(simulate_without_showers("TFCSBinnedShower", energy), FCSFatal);
  EXPECT_DOUBLE_EQ(energy, 0.0);
}

// Without a loaded model, TFCSBinnedShowerONNX fails instead of dereferencing
// the missing simulator
TEST(BinnedShowerTests, SimulateWithoutModel)
{
  double energy = -1;
  EXPECT_EQ(simulate_without_showers("TFCSBinnedShowerONNX", energy), FCSFatal);
  EXPECT_DOUBLE_EQ(energy, 0.0);
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
