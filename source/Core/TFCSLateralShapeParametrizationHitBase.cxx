// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#include "FastCaloSim/Core/TFCSLateralShapeParametrizationHitBase.h"

#include "FastCaloSim/Core/TFCSExtrapolationState.h"
#include "FastCaloSim/Core/TFCSSimulationState.h"
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
    const TFCSExtrapolationState* extrapol) const -> float
{
  const int nhits = get_number_of_hits(simulstate, truth, extrapol);
  const int sample = calosample();
  if (nhits <= 0 || sample < 0)
    return -1.;
  else
    return simulstate.E(sample) / nhits;
}

auto TFCSLateralShapeParametrizationHitBase::getMinWeight() const -> float
{
  return -1.;
}

auto TFCSLateralShapeParametrizationHitBase::getMaxWeight() const -> float
{
  return -1.;
}

auto TFCSLateralShapeParametrizationHitBase::simulate_hit(
    Hit& hit,
    TFCSSimulationState& /*simulstate*/,
    const TFCSTruthState* /*truth*/,
    const TFCSExtrapolationState* extrapol) -> FCSReturnCode
{
  int cs = calosample();
  hit.set_eta_x(0.5
                * (extrapol->eta(cs, Cell::SubPos::ENT)
                   + extrapol->eta(cs, Cell::SubPos::EXT)));
  hit.set_phi_y(0.5
                * (extrapol->phi(cs, Cell::SubPos::ENT)
                   + extrapol->phi(cs, Cell::SubPos::EXT)));

  return FCSSuccess;
}
