/*
  Copyright (C) 2002-2024 CERN for the benefit of the ATLAS collaboration
*/

#include <utility>

#include "FastCaloSim/Core/TFCSHistoLateralShapeWeight.h"

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGaussZiggurat.h"
#include "FastCaloSim/Core/TFCSSimulationState.h"
#include "TH1.h"
#include "TMath.h"
#include "TVector2.h"

//=============================================
//======= TFCSHistoLateralShapeWeight =========
//=============================================

TFCSHistoLateralShapeWeight::TFCSHistoLateralShapeWeight(const char* name,
                                                         const char* title)
    : TFCSLateralShapeParametrizationHitBase(name, title)
{
}

TFCSHistoLateralShapeWeight::~TFCSHistoLateralShapeWeight()
{
  if (m_hist)
    delete m_hist;
}

float TFCSHistoLateralShapeWeight::getMinWeight() const
{
  return m_minWeight;
}

float TFCSHistoLateralShapeWeight::getMaxWeight() const
{
  return m_maxWeight;
}

FCSReturnCode TFCSHistoLateralShapeWeight::simulate_hit(
    Hit& hit,
    TFCSSimulationState& simulstate,
    const TFCSTruthState* /*truth*/,
    const TFCSExtrapolationState* /*extrapol*/)
{
  if (!simulstate.randomEngine()) {
    return FCSFatal;
  }

  const double center_eta = hit.center_eta();
  const double center_phi = hit.center_phi();
  const double center_r = hit.center_r();
  const double center_z = hit.center_z();

  const float dist000 = TMath::Sqrt(center_r * center_r + center_z * center_z);
  const float eta_jakobi = TMath::Abs(2.0 * TMath::Exp(-center_eta)
                                      / (1.0 + TMath::Exp(-2 * center_eta)));

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
    weight =
        CLHEP::RandGaussZiggurat::shoot(simulstate.randomEngine(), weight, RMS);
  }
  hit.E() *= weight;

  ATH_MSG_DEBUG("HIT: E=" << hit.E() << " dR_mm=" << delta_r_mm
                          << " weight=" << weight);
  return FCSSuccess;
}

bool TFCSHistoLateralShapeWeight::Initialize(TH1* hist)
{
  if (!hist)
    return false;
  if (m_hist)
    delete m_hist;
  m_hist = (TH1*)hist->Clone(TString("TFCSHistoLateralShapeWeight_")
                             + hist->GetName());
  m_hist->SetDirectory(nullptr);

  return true;
}

void TFCSHistoLateralShapeWeight::Print(Option_t* option) const
{
  TString opt(option);
  bool shortprint = opt.Index("short") >= 0;
  bool longprint = msgLvl(MSG::DEBUG) || (msgLvl(MSG::INFO) && !shortprint);
  TString optprint = opt;
  optprint.ReplaceAll("short", "");
  TFCSLateralShapeParametrizationHitBase::Print(option);

  if (longprint) {
    if (m_hist)
      ATH_MSG_INFO(optprint
                   << "  Histogram: " << m_hist->GetNbinsX() << " bins ["
                   << std::as_const(m_hist)->GetXaxis()->GetXmin() << ","
                   << std::as_const(m_hist)->GetXaxis()->GetXmax() << "]"
                   << " min weight: " << m_minWeight
                   << " max weight: " << m_maxWeight);
    else
      ATH_MSG_INFO(optprint << "  no Histogram");
  }
}
