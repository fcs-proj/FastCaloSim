/*
  Copyright (C) 2002-2019 CERN for the benefit of the ATLAS collaboration
*/

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGaussZiggurat.h"

#include "FastCaloSim/TFCSHistoLateralShapeWeightHitAndMiss.h"
#include "FastCaloSim/TFCSSimulationState.h"

#include "TH1.h"
#include "TVector2.h"
#include "TMath.h"

//=============================================
//======= TFCSHistoLateralShapeWeightHitAndMiss =========
//=============================================

TFCSHistoLateralShapeWeightHitAndMiss::TFCSHistoLateralShapeWeightHitAndMiss(
    const char *name, const char *title)
    : TFCSHistoLateralShapeWeight(name, title) {}

TFCSHistoLateralShapeWeightHitAndMiss::
    ~TFCSHistoLateralShapeWeightHitAndMiss() {}

FCSReturnCode TFCSHistoLateralShapeWeightHitAndMiss::simulate_hit(
    Hit &hit, TFCSSimulationState &simulstate, const TFCSTruthState * /*truth*/,
    const TFCSExtrapolationState * /*extrapol*/) {
  if (!simulstate.randomEngine()) {
    return FCSFatal;
  }

  const double center_eta = hit.center_eta();
  const double center_phi = hit.center_phi();
  const double center_r = hit.center_r();
  const double center_z = hit.center_z();

  const float dist000 = TMath::Sqrt(center_r * center_r + center_z * center_z);
  const float eta_jakobi = TMath::Abs(2.0 * TMath::Exp(-center_eta) /
                                      (1.0 + TMath::Exp(-2 * center_eta)));

  const float delta_eta = hit.eta() - center_eta;
  const float delta_phi = hit.phi() - center_phi;
  const float delta_eta_mm = delta_eta * eta_jakobi * dist000;
  const float delta_phi_mm = delta_phi * center_r;
  const float delta_r_mm =
      TMath::Sqrt(delta_eta_mm * delta_eta_mm + delta_phi_mm * delta_phi_mm);

  // TODO: delta_r_mm should perhaps be cached in hit

  Int_t bin = m_hist->FindBin(delta_r_mm);
  if (bin < 1)
    bin = 1;
  if (bin > m_hist->GetNbinsX())
    bin = m_hist->GetNbinsX();
  float meanweight = m_hist->GetBinContent(bin);
  float weight = meanweight;
  float RMS = m_hist->GetBinError(bin);
  if (RMS > 0) {
    weight = CLHEP::RandGaussZiggurat::shoot(simulstate.randomEngine(),
                                             meanweight, RMS);
  }

  /* -------------------------------------------------------------------
   * Weight is used to scale hit energy.
   *
   * if (meanweight > m_minWeight and meanweight < m_maxWeight)
   * 	Hit is accecpted with probability of m_minWeight/meanweight.
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
  if (meanweight >= m_maxWeight) {
    weight = m_maxWeight;
  } else if (meanweight > m_minWeight) {
    float prob = m_minWeight / meanweight;
    float rnd = CLHEP::RandFlat::shoot(simulstate.randomEngine());
    if (rnd >= prob)
      weight = 0.;
  }

  hit.E() *= weight;
  ATH_MSG_DEBUG("HIT: E=" << hit.E() << " dR_mm=" << delta_r_mm
                          << " meanweight=" << meanweight
                          << " weight=" << weight);

  return FCSSuccess;
}
