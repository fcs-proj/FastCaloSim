// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#include <algorithm>
#include <cmath>

#include "FastCaloSim/Core/TFCSHistoLateralShapeGausLogWeightHitAndMiss.h"

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGaussZiggurat.h"
#include "FastCaloSim/Core/GaussZiggurat.h"
#include "FastCaloSim/Core/TFCSSimulationState.h"
#include "FastCaloSim/Definitions/Phi.h"
#include "TH1.h"
#include "TVector2.h"

//=============================================
//======= TFCSHistoLateralShapeGausLogWeightHitAndMiss =========
//=============================================

TFCSHistoLateralShapeGausLogWeightHitAndMiss::
    TFCSHistoLateralShapeGausLogWeightHitAndMiss(const char* name,
                                                 const char* title)
    : TFCSHistoLateralShapeWeight(name, title)
{
}

TFCSHistoLateralShapeGausLogWeightHitAndMiss::
    ~TFCSHistoLateralShapeGausLogWeightHitAndMiss() = default;

auto TFCSHistoLateralShapeGausLogWeightHitAndMiss::simulate_hit(
    Hit& hit,
    TFCSSimulationState& simulstate,
    const TFCSTruthState* /*truth*/,
    const TFCSExtrapolationState* /*extrapol*/) -> FCSReturnCode
{
  if (simulstate.randomEngine() == nullptr) {
    return FCSFatal;
  }

  const double center_eta = hit.center_eta();
  const double center_phi = hit.center_phi();
  const double center_r = hit.center_r();
  const double center_z = hit.center_z();

  const double dist000 = FastCaloSim::dist000(center_r, center_z);
  const double eta_jakobi = FastCaloSim::eta_jakobi(center_eta);

  const double delta_eta = hit.eta() - center_eta;
  const double delta_phi = FastCaloSim::delta_phi(hit.phi(), center_phi);
  const double delta_eta_mm = delta_eta * eta_jakobi * dist000;
  const double delta_phi_mm = delta_phi * center_r;
  const double delta_r_mm =
      std::sqrt(delta_eta_mm * delta_eta_mm + delta_phi_mm * delta_phi_mm);

  // TODO: delta_r_mm should perhaps be cached in hit

  Int_t bin = m_hist->FindBin(delta_r_mm);
  bin = std::max(bin, 1);
  bin = std::min(bin, m_hist->GetNbinsX());
  double const meanweight = m_hist->GetBinContent(bin);
  double weight = meanweight;
  double const RMS = m_hist->GetBinError(bin);
  if (RMS > 0) {
    double const logweight = FastCaloSim::gauss_ziggurat(
        simulstate.randomEngine(), -0.5 * RMS * RMS, RMS);
    weight *= std::exp(logweight);
  }
  if (meanweight <= 1) {
    // if meanweight<=1, give lower energy to hit.
    // TFCSLateralShapeParametrizationHitChain needs to be able to generate more
    // hits in this case
    hit.set_E(weight * hit.E());
  } else {
    // if meanweight>1, accept only 1/meanweight events, but give them a higher
    // energy increased by weight. This leads to larger fluctuations, while
    // keeping the shape unchanged.
    double const prob = 1.0 / meanweight;
    double const rnd = CLHEP::RandFlat::shoot(simulstate.randomEngine());
    if (rnd < prob) {
      hit.set_E(weight * hit.E());
    } else {
      hit.set_E(0);
    }
  }

  FCS_MSG_DEBUG("HIT: E=" << hit.E() << " dR_mm=" << delta_r_mm
                          << " meanweight=" << meanweight
                          << " weight=" << weight);
  return FCSSuccess;
}
