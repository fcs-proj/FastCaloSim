/*
  Copyright (C) 2002-2018 CERN for the benefit of the ATLAS collaboration
*/

#include "FastCaloSim/TFCSLateralShapeParametrizationHitBase.h"

#include "FastCaloSim/FastCaloSim_CaloCell_ID.h"
#include "FastCaloSim/TFCSExtrapolationState.h"
#include "FastCaloSim/TFCSSimulationState.h"

//=============================================
//======= TFCSLateralShapeParametrization =========
//=============================================

TFCSLateralShapeParametrizationHitBase::TFCSLateralShapeParametrizationHitBase(
    const char* name, const char* title)
    : TFCSLateralShapeParametrization(name, title)
{
}

double TFCSLateralShapeParametrizationHitBase::get_sigma2_fluctuation(
    TFCSSimulationState& /*simulstate*/,
    const TFCSTruthState* /*truth*/,
    const TFCSExtrapolationState* /*extrapol*/) const
{
  return -1;
}

int TFCSLateralShapeParametrizationHitBase::get_number_of_hits(
    TFCSSimulationState& /*simulstate*/,
    const TFCSTruthState* /*truth*/,
    const TFCSExtrapolationState* /*extrapol*/) const
{
  return -1;
}

float TFCSLateralShapeParametrizationHitBase::get_E_hit(
    TFCSSimulationState& simulstate,
    const TFCSTruthState* truth,
    const TFCSExtrapolationState* extrapol) const
{
  const int nhits = get_number_of_hits(simulstate, truth, extrapol);
  const int sample = calosample();
  if (nhits <= 0 || sample < 0)
    return -1.;
  else
    return simulstate.E(sample) / nhits;
}

float TFCSLateralShapeParametrizationHitBase::getMinWeight() const
{
  return -1.;
}

float TFCSLateralShapeParametrizationHitBase::getMaxWeight() const
{
  return -1.;
}

FCSReturnCode TFCSLateralShapeParametrizationHitBase::simulate_hit(
    Hit& hit,
    TFCSSimulationState& /*simulstate*/,
    const TFCSTruthState* /*truth*/,
    const TFCSExtrapolationState* extrapol)
{
  int cs = calosample();
  hit.eta() = 0.5
      * (extrapol->eta(cs, CaloSubPos::SUBPOS_ENT)
         + extrapol->eta(cs, CaloSubPos::SUBPOS_EXT));
  hit.phi() = 0.5
      * (extrapol->phi(cs, CaloSubPos::SUBPOS_ENT)
         + extrapol->phi(cs, CaloSubPos::SUBPOS_EXT));

  return FCSSuccess;
}
