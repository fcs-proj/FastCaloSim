// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#include <cstdint>
#include <memory>

#include <CLHEP/Random/MixMaxRng.h>
#include <TClass.h>
#include <gtest/gtest.h>

#include "FastCaloSim/Core/TFCSLateralShapeParametrization.h"
#include "FastCaloSim/Core/TFCSLateralShapeParametrizationFluctChain.h"
#include "FastCaloSim/Core/TFCSParametrizationBase.h"
#include "FastCaloSim/Core/TFCSSimulationState.h"
#include "FastCaloSim/Core/TFCSTruthState.h"
#include "FastCaloSim/Geometry/CaloGeo.h"
#include "FastCaloSim/Geometry/Cell.h"

namespace
{
// Created through ROOT, as when read from a parametrization file. These
// classes are not exported, so the test owns them through the exported base
// class.
auto make(const char* class_name, std::size_t layer)
    -> std::unique_ptr<TFCSParametrizationBase>
{
  std::unique_ptr<TFCSParametrizationBase> param(
      static_cast<TFCSParametrizationBase*>(
          TClass::GetClass(class_name)->New()));
  static_cast<TFCSLateralShapeParametrization*>(param.get())
      ->set_calosample(layer);
  return param;
}

// A geometry whose only cell is far from every hit: TFCSHitCellMapping
// rejects each hit, keeps it with zero energy and reports success
class FarCellGeo : public CaloGeo
{
public:
  auto get_cell(std::size_t /*layer*/, const Position& /*pos*/) const
      -> Cell override
  {
    const Position center {
        .m_x = 0, .m_y = 0, .m_z = 0, .m_eta = 3.0, .m_phi = 0, .m_r = 0};
    return {1,
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
            0.1,
            0.1,
            0};
  }
};
}  // namespace

// The chain ends when its hits never deposit energy, instead of simulating
// hits forever
TEST(FluctChainTest, EndsWhenNoHitDepositsEnergy)
{
  constexpr std::size_t layer = 0;
  FarCellGeo geo;
  const auto chain = make("TFCSLateralShapeParametrizationFluctChain", layer);
  const auto number_of_hits =
      make("TFCSLateralShapeParametrizationHitNumberFromE", layer);
  const auto mapping = make("TFCSHitCellMapping", layer);
  mapping->set_geometry(&geo);

  auto& fluct_chain =
      *static_cast<TFCSLateralShapeParametrizationFluctChain*>(chain.get());
  fluct_chain.set_number_of_hits_simul(
      static_cast<TFCSLateralShapeParametrizationHitBase*>(
          number_of_hits.get()));
  fluct_chain.push_back(
      static_cast<TFCSLateralShapeParametrizationHitBase*>(mapping.get()));

  CLHEP::MixMaxRng engine(42);
  TFCSSimulationState simulstate(&engine);
  simulstate.set_E(layer, 1000.0);
  simulstate.set_E(1000.0);
  TFCSTruthState truth;
  truth.set_pdgid(22);

  EXPECT_EQ(chain->simulate(simulstate, &truth, nullptr), FCSSuccess);
  EXPECT_TRUE(simulstate.cells().empty());
}
