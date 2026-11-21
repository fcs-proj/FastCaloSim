// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#include <algorithm>
#include <cmath>

#include "FastCaloSim/Core/TFCSHistoLateralShapeWeightHitAndMiss.h"

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGaussZiggurat.h"
#include "FastCaloSim/Core/GaussZiggurat.h"
#include "FastCaloSim/Core/TFCSSimulationState.h"
#include "TH1.h"
#include "TVector2.h"

//=============================================
//======= TFCSHistoLateralShapeWeightHitAndMiss =========
//=============================================

TFCSHistoLateralShapeWeightHitAndMiss::TFCSHistoLateralShapeWeightHitAndMiss(
    const char* name, const char* title)
    : TFCSHistoLateralShapeWeight(name, title)
{
}

TFCSHistoLateralShapeWeightHitAndMiss::
    ~TFCSHistoLateralShapeWeightHitAndMiss() = default;

auto TFCSHistoLateralShapeWeightHitAndMiss::simulate_hit(
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
  const double delta_phi = hit.phi() - center_phi;
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
    weight =
        FastCaloSim::gauss_ziggurat(simulstate.randomEngine(), meanweight, RMS);
  }

  /* -------------------------------------------------------------------
   * Weight is used to scale hit energy.
   *
   * if (meanweight > m_minWeight and meanweight < m_maxWeight)
   * 	Hit is accepted with probability of m_minWeight/meanweight.
   * 	If not accepted, weight is set to zero (this is
   *    equivalent to not accept the hit).
   *
   * if (meanweight<=m_minWeight)
   * 	Give lower energy to hit with probability 1.
   * 	TFCSLateralShapeParametrizationHitChain needs to be able
   * 	to generate more hits in this case
   *
   * if (meanweight >= m_maxWeight)
   * 	meanweight is above upper threshold. It is set to m_maxWeight.
   * -------------------------------------------------------------------
   */
  if (meanweight >= static_cast<double>(m_maxWeight)) {
    weight = static_cast<double>(m_maxWeight);
  } else if (meanweight > static_cast<double>(m_minWeight)) {
    double const prob = static_cast<double>(m_minWeight) / meanweight;
    double const rnd = CLHEP::RandFlat::shoot(simulstate.randomEngine());
    if (rnd >= prob) {
      weight = 0.;
    }
  }

  hit.set_E(weight * hit.E());
  FCS_MSG_DEBUG("HIT: E=" << hit.E() << " dR_mm=" << delta_r_mm
                          << " meanweight=" << meanweight
                          << " weight=" << weight);

  return FCSSuccess;
}
