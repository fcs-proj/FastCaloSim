// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#include "FastCaloSim/Core/TFCSLateralShapeParametrizationHitBase.h"

#include "FastCaloSim/Core/TFCSExtrapolationState.h"
#include "FastCaloSim/Core/TFCSSimulationState.h"
#include "FastCaloSim/Definitions/Phi.h"
#include "FastCaloSim/Geometry/Cell.h"

//=============================================
//======= TFCSLateralShapeParametrization =========
//=============================================

TFCSLateralShapeParametrizationHitBase::TFCSLateralShapeParametrizationHitBase(
    const char* name, const char* title)
    : TFCSLateralShapeParametrization(name, title)
{
}

auto TFCSLateralShapeParametrizationHitBase::get_sigma2_fluctuation(
    TFCSSimulationState& /*simulstate*/,
    const TFCSTruthState* /*truth*/,
    const TFCSExtrapolationState* /*extrapol*/) const -> double
{
  return -1;
}

auto TFCSLateralShapeParametrizationHitBase::get_number_of_hits(
    TFCSSimulationState& /*simulstate*/,
    const TFCSTruthState* /*truth*/,
    const TFCSExtrapolationState* /*extrapol*/) const -> int
{
  return -1;
}

auto TFCSLateralShapeParametrizationHitBase::get_E_hit(
    TFCSSimulationState& simulstate,
    const TFCSTruthState* truth,
    const TFCSExtrapolationState* extrapol) const -> double
{
  const int nhits = get_number_of_hits(simulstate, truth, extrapol);
  if (nhits <= 0 || !has_calosample()) {
    return -1.;
  }
  return simulstate.E(calosample()) / nhits;
}

auto TFCSLateralShapeParametrizationHitBase::getMinWeight() const -> double
{
  return -1.;
}

auto TFCSLateralShapeParametrizationHitBase::getMaxWeight() const -> double
{
  return -1.;
}

auto TFCSLateralShapeParametrizationHitBase::simulate_hit(
    Hit& hit,
    TFCSSimulationState& /*simulstate*/,
    const TFCSTruthState* /*truth*/,
    const TFCSExtrapolationState* extrapol) -> FCSReturnCode
{
  std::size_t const cs = calosample();
  hit.set_eta_x(0.5
                * (extrapol->eta(cs, FastCaloSim::SubPos::ENT)
                   + extrapol->eta(cs, FastCaloSim::SubPos::EXT)));
  hit.set_phi_y(
      FastCaloSim::interpolate_phi(extrapol->phi(cs, FastCaloSim::SubPos::ENT),
                                   extrapol->phi(cs, FastCaloSim::SubPos::EXT),
                                   0.5));

  return FCSSuccess;
}
