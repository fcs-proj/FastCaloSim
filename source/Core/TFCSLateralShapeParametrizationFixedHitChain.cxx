// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#include "FastCaloSim/Core/TFCSLateralShapeParametrizationFixedHitChain.h"

#include "CLHEP/Random/RandGaussZiggurat.h"
#include "FastCaloSim/Core/TFCSSimulationState.h"
#include "TMath.h"

//=============================================
//======= TFCSLateralShapeParametrizationFixedHitChain =========
//=============================================

TFCSLateralShapeParametrizationFixedHitChain::
    TFCSLateralShapeParametrizationFixedHitChain(const char* name,
                                                 const char* title)
    : TFCSLateralShapeParametrizationHitChain(name, title)
{
}

TFCSLateralShapeParametrizationFixedHitChain::
    TFCSLateralShapeParametrizationFixedHitChain(
        TFCSLateralShapeParametrizationHitBase* hitsim)
    : TFCSLateralShapeParametrizationHitChain(hitsim)
{
}

bool TFCSLateralShapeParametrizationFixedHitChain::check_all_hits_simulated(
    TFCSLateralShapeParametrizationHitBase::Hit& hit,
    TFCSSimulationState& simulstate,
    const TFCSTruthState* truth,
    const TFCSExtrapolationState* extrapol,
    bool success) const
{
  (void)success;  // unused parameter
  (void)truth;  // unused parameter
  (void)extrapol;  // unused parameter

  if (!simulstate.hasAuxInfo("FCSHitChainNHits"_FCShash)) {
    FCS_MSG_ERROR("FCSHitChainNHits aux-info not set - aborting simulation!");
    return true;
  }

  int nhit_signed = simulstate.getAuxInfo<int>("FCSHitChainNHits"_FCShash);

  if (nhit_signed < 0) {
    FCS_MSG_ERROR(
        "The number of hits to simulate is not defined - aborting simulation!");
    return true;
  }

  unsigned int nhit = static_cast<unsigned int>(nhit_signed);

  long unsigned int current_hit = hit.idx();

  // Take care of layer energy updating.
  // Needed for the phi modulation!
  bool done = current_hit >= nhit;
  if (!done) {
    double sumEhit = 0;
    if (simulstate.hasAuxInfo("FCSFixedHitChainEnergySum"_FCShash)
        && hit.idx() > 0)
    {
      sumEhit =
          simulstate.getAuxInfo<double>("FCSFixedHitChainEnergySum"_FCShash);
    }

    sumEhit += hit.E();
    simulstate.setAuxInfo<double>("FCSFixedHitChainEnergySum"_FCShash, sumEhit);
  } else {
    if (!simulstate.hasAuxInfo("FCSFixedHitChainEnergySum"_FCShash)) {
      simulstate.setAuxInfo<double>("FCSFixedHitChainEnergySum"_FCShash, 0.);
    }

    double old_energy = simulstate.E(calosample());
    double new_energy =
        simulstate.getAuxInfo<double>("FCSFixedHitChainEnergySum"_FCShash);

    double energy_difference = new_energy - old_energy;

    simulstate.set_E(simulstate.E() + energy_difference);
    simulstate.set_E(calosample(), new_energy);

    if (simulstate.E() > std::numeric_limits<double>::epsilon()) {
      simulstate.set_Efrac(calosample(), new_energy / simulstate.E());
    }
    simulstate.setAuxInfo<double>("FCSFixedHitChainEnergySum"_FCShash, 0.);
  }

  return done;
}
