/*
  Copyright (C) 2002-2019 CERN for the benefit of the ATLAS collaboration
*/

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGaussZiggurat.h"

#include "FastCaloSim/TFCSHistoLateralShapeGausLogWeightHitAndMiss.h"
#include "FastCaloSim/TFCSSimulationState.h"

#include "TH1.h"
#include "TVector2.h"
#include "TMath.h"

//=============================================
//======= TFCSHistoLateralShapeGausLogWeightHitAndMiss =========
//=============================================

TFCSHistoLateralShapeGausLogWeightHitAndMiss::
    TFCSHistoLateralShapeGausLogWeightHitAndMiss(const char *name,
                                                 const char *title)
    : TFCSHistoLateralShapeWeight(name, title) {}

TFCSHistoLateralShapeGausLogWeightHitAndMiss::
    ~TFCSHistoLateralShapeGausLogWeightHitAndMiss() {}

FCSReturnCode TFCSHistoLateralShapeGausLogWeightHitAndMiss::simulate_hit(
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
    float logweight = CLHEP::RandGaussZiggurat::shoot(simulstate.randomEngine(),
                                                      -0.5 * RMS * RMS, RMS);
    weight *= TMath::Exp(logweight);
  }
  if (meanweight <= 1) {
    // if meanweight<=1, give lower energy to hit.
    // TFCSLateralShapeParametrizationHitChain needs to be able to generate more
    // hits in this case
    hit.E() *= weight;
  } else {
    // if meanweight>1, accept only 1/meanweight events, but give them a higher
    // energy increased by weight. This leads to larger fluctuations, while
    // keeping the shape unchanged.
    float prob = 1.0 / meanweight;
    float rnd = CLHEP::RandFlat::shoot(simulstate.randomEngine());
    if (rnd < prob)
      hit.E() *= weight;
    else
      hit.E() = 0;
  }

  ATH_MSG_DEBUG("HIT: E=" << hit.E() << " dR_mm=" << delta_r_mm
                          << " meanweight=" << meanweight
                          << " weight=" << weight);
  return FCSSuccess;
}
