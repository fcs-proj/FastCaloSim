// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#include <cstddef>
#include <limits>

#include "FastCaloSim/Core/TFCSLateralShapeParametrizationFixedHitChain.h"

#include "CLHEP/Random/RandGaussZiggurat.h"
#include "FastCaloSim/Core/TFCSSimulationState.h"

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

auto TFCSLateralShapeParametrizationFixedHitChain::check_all_hits_simulated(
    TFCSLateralShapeParametrizationHitBase::Hit& hit,
    TFCSSimulationState& simulstate,
    const TFCSTruthState* truth,
    const TFCSExtrapolationState* extrapol,
    bool success,
    double sum_E_hits) const -> bool
{
  (void)truth;  // unused parameter
  (void)extrapol;  // unused parameter

  if (!success) {
    // The hit simulation failed => Nothing changed => Still not done
    return false;
  }

  if (!simulstate.hasAuxInfo("FCSHitChainNHits"_FCShash)) {
    FCS_MSG_ERROR("FCSHitChainNHits aux-info not set - aborting simulation!");
    return true;
  }

  int const nhit_signed =
      simulstate.getAuxInfo<int>("FCSHitChainNHits"_FCShash);

  if (nhit_signed < 0) {
    FCS_MSG_ERROR(
        "The number of hits to simulate is not defined - aborting simulation!");
    return true;
  }

  if (hit.idx() < static_cast<std::size_t>(nhit_signed)) {
    return false;
  }

  // All hits are simulated: the layer energy becomes the energy of these hits
  // (needed for the phi modulation)
  double const old_energy = simulstate.E(calosample());
  double const energy_difference = sum_E_hits - old_energy;

  simulstate.set_E(simulstate.E() + energy_difference);
  simulstate.set_E(calosample(), sum_E_hits);

  if (simulstate.E() > std::numeric_limits<double>::epsilon()) {
    simulstate.set_Efrac(calosample(), sum_E_hits / simulstate.E());
  }
  return true;
}
