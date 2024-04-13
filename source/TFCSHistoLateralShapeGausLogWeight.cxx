/*
  Copyright (C) 2002-2019 CERN for the benefit of the ATLAS collaboration
*/

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGaussZiggurat.h"

#include "FastCaloSim/TFCSHistoLateralShapeGausLogWeight.h"
#include "FastCaloSim/TFCSSimulationState.h"

#include "TH1.h"
#include "TVector2.h"
#include "TMath.h"

//=============================================
//======= TFCSHistoLateralShapeGausLogWeight =========
//=============================================

TFCSHistoLateralShapeGausLogWeight::TFCSHistoLateralShapeGausLogWeight(
    const char *name, const char *title)
    : TFCSHistoLateralShapeWeight(name, title) {}

TFCSHistoLateralShapeGausLogWeight::~TFCSHistoLateralShapeGausLogWeight() {
#ifdef USE_GPU
  delete m_LdFH;
#endif
}

FCSReturnCode TFCSHistoLateralShapeGausLogWeight::simulate_hit(
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
  float weight = m_hist->GetBinContent(bin);
  float RMS = m_hist->GetBinError(bin);
  if (RMS > 0) {
    float logweight = CLHEP::RandGaussZiggurat::shoot(simulstate.randomEngine(),
                                                      -0.5 * RMS * RMS, RMS);
    weight *= TMath::Exp(logweight);
  }
  hit.E() *= weight;

  ATH_MSG_DEBUG("HIT: E=" << hit.E() << " dR_mm=" << delta_r_mm
                          << " weight=" << weight);
  return FCSSuccess;
}

#ifdef USE_GPU
// for FCS-GPU
// make the histograms with GPU EDM
void TFCSHistoLateralShapeGausLogWeight::LoadHist() {

  if (m_LdFH) {
    return;
  }
  m_LdFH = new LoadGpuHist();
  FH1D fh1d;

  fh1d.nbins = m_hist->GetNbinsX();
  fh1d.h_contents = (float *)malloc(fh1d.nbins * sizeof(float));
  fh1d.h_errors = (float *)malloc(fh1d.nbins * sizeof(float));
  fh1d.h_borders = (float *)malloc((fh1d.nbins + 1) * sizeof(float));

  for (int ibin = 0; ibin < fh1d.nbins; ++ibin) {
    fh1d.h_contents[ibin] = m_hist->GetBinContent(ibin + 1);
    fh1d.h_errors[ibin] = m_hist->GetBinError(ibin + 1);
    fh1d.h_borders[ibin] = m_hist->GetXaxis()->GetBinLowEdge(ibin + 1);
  }
  fh1d.h_borders[fh1d.nbins] = m_hist->GetXaxis()->GetBinUpEdge(fh1d.nbins);

  m_LdFH->set_hf1d(&fh1d);
  m_LdFH->LD1D();

  return;
}
#endif
