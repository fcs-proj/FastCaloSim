// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#include <algorithm>
#include <cmath>

#include "FastCaloSim/Core/TFCSHistoLateralShapeGausLogWeight.h"

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGaussZiggurat.h"
#include "FastCaloSim/Core/GaussZiggurat.h"
#include "FastCaloSim/Core/TFCSSimulationState.h"
#include "TH1.h"
#include "TVector2.h"

//=============================================
//======= TFCSHistoLateralShapeGausLogWeight =========
//=============================================

TFCSHistoLateralShapeGausLogWeight::TFCSHistoLateralShapeGausLogWeight(
    const char* name, const char* title)
    : TFCSHistoLateralShapeWeight(name, title)
{
}

TFCSHistoLateralShapeGausLogWeight::~TFCSHistoLateralShapeGausLogWeight() =
    default;

auto TFCSHistoLateralShapeGausLogWeight::simulate_hit(
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
  double weight = m_hist->GetBinContent(bin);
  double const RMS = m_hist->GetBinError(bin);
  if (RMS > 0) {
    double const logweight = FastCaloSim::gauss_ziggurat(
        simulstate.randomEngine(), -0.5 * RMS * RMS, RMS);
    weight *= std::exp(logweight);
  }
  hit.set_E(weight * hit.E());

  FCS_MSG_DEBUG("HIT: E=" << hit.E() << " dR_mm=" << delta_r_mm
                          << " weight=" << weight);
  return FCSSuccess;
}
