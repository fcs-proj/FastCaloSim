// Copyright (c) 2024 CERN for the benefit of the FastCaloSim project

#include <cmath>

#include "FastCaloSim/Core/TFCSHistoLateralShapeParametrization.h"

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandPoisson.h"
#include "FastCaloSim/Core/TFCSExtrapolationState.h"
#include "FastCaloSim/Core/TFCSSimulationState.h"
#include "FastCaloSim/Definitions/ParticleData.h"
#include "TFile.h"
#include "TH2.h"
#include "TMath.h"

//=============================================
//======= TFCSHistoLateralShapeParametrization =========
//=============================================

TFCSHistoLateralShapeParametrization::TFCSHistoLateralShapeParametrization(
    const char* name, const char* title)
    : TFCSLateralShapeParametrizationHitBase(name, title)
    , m_nhits(0)
    , m_r_offset(0)
    , m_r_scale(1.0)
{
  TFCSHistoLateralShapeParametrization::reset_phi_symmetric();
}

TFCSHistoLateralShapeParametrization::~TFCSHistoLateralShapeParametrization() {}

void TFCSHistoLateralShapeParametrization::set_geometry(CaloGeo* geo)
{
  TFCSLateralShapeParametrizationHitBase::set_geometry(geo);
  if (!m_hist.get_HistoContents().empty()) {
    int first_fix_bin = -1;
    for (int i = (int)(m_hist.get_HistoContents().size() - 1); i >= 0; --i) {
      if (std::isnan(m_hist.get_HistoContents()[i])) {
        FCS_MSG_DEBUG("nan in histo content for "
                      << GetTitle() << ", bin[" << i
                      << "]=" << m_hist.get_HistoContents()[i] << " -> 1");
        m_hist.get_HistoContents()[i] = 1;
        first_fix_bin = i;
      }
    }
    if (first_fix_bin < 0)
      return;

    if (first_fix_bin == 0) {
      FCS_MSG_WARNING("nan in histo content for "
                      << GetTitle()
                      << " for all bins. Fixed to probability 1 causing hits "
                         "to be deposited in the shower center");
    } else {
      int last_fix_bin = -1;
      for (size_t i = 0; i < m_hist.get_HistoContents().size(); ++i) {
        if (std::isnan(m_hist.get_HistoContents()[i])) {
          FCS_MSG_DEBUG("nan in histo content for "
                        << GetTitle() << ", bin[" << i
                        << "]=" << m_hist.get_HistoContents()[i] << " -> 0");
          m_hist.get_HistoContents()[i] = 0;
          last_fix_bin = i;
        }
      }
      FCS_MSG_WARNING("nan in histo content for "
                      << GetTitle() << ". Fixed up to bin " << last_fix_bin
                      << " with probability 0 and beyond bin " << first_fix_bin
                      << " with probability 1.");
    }
  }
}

double TFCSHistoLateralShapeParametrization::get_sigma2_fluctuation(
    TFCSSimulationState& /*simulstate*/,
    const TFCSTruthState* /*truth*/,
    const TFCSExtrapolationState* /*extrapol*/) const
{
  // Limit to factor 1000 fluctuations
  if (m_nhits < 0.001)
    return 1000;
  return 1.0 / m_nhits;
}

int TFCSHistoLateralShapeParametrization::get_number_of_hits(
    TFCSSimulationState& simulstate,
    const TFCSTruthState* /*truth*/,
    const TFCSExtrapolationState* /*extrapol*/) const
{
  if (!simulstate.randomEngine()) {
    return -1;
  }

  return CLHEP::RandPoisson::shoot(simulstate.randomEngine(), m_nhits);
}

void TFCSHistoLateralShapeParametrization::set_number_of_hits(float nhits)
{
  m_nhits = nhits;
}

FCSReturnCode TFCSHistoLateralShapeParametrization::simulate_hit(
    Hit& hit,
    TFCSSimulationState& simulstate,
    const TFCSTruthState* truth,
    const TFCSExtrapolationState* /*extrapol*/)
{
  if (!simulstate.randomEngine()) {
    return FCSFatal;
  }

  const int pdgId = truth->pdgid();
  const double charge = ParticleData::charge(pdgId);

  const int cs = calosample();
  const double center_eta = hit.center_eta();
  const double center_phi = hit.center_phi();
  const double center_r = hit.center_r();
  const double center_z = hit.center_z();

  if (TMath::IsNaN(center_r) or TMath::IsNaN(center_z)
      or TMath::IsNaN(center_eta) or TMath::IsNaN(center_phi))
  {  // Check if extrapolation fails
    return FCSFatal;
  }

  float alpha, r, rnd1, rnd2;
  rnd1 = CLHEP::RandFlat::shoot(simulstate.randomEngine());
  rnd2 = CLHEP::RandFlat::shoot(simulstate.randomEngine());
  if (is_phi_symmetric()) {
    if (rnd2 >= 0.5) {  // Fill negative phi half of shape
      rnd2 -= 0.5;
      rnd2 *= 2;
      m_hist.rnd_to_fct(alpha, r, rnd1, rnd2);
      alpha = -alpha;
    } else {  // Fill positive phi half of shape
      rnd2 *= 2;
      m_hist.rnd_to_fct(alpha, r, rnd1, rnd2);
    }
  } else {
    m_hist.rnd_to_fct(alpha, r, rnd1, rnd2);
  }
  if (TMath::IsNaN(alpha) || TMath::IsNaN(r)) {
    FCS_MSG_ERROR("  Histogram: "
                  << m_hist.get_HistoBordersx().size() - 1 << "*"
                  << m_hist.get_HistoBordersy().size() - 1
                  << " bins, #hits=" << m_nhits << " alpha=" << alpha
                  << " r=" << r << " rnd1=" << rnd1 << " rnd2=" << rnd2);
    alpha = 0;
    r = 0.001;

    FCS_MSG_ERROR("  This error could probably be retried");
    return FCSFatal;
  }

  r *= m_r_scale;
  r += m_r_offset;
  if (r < 0)
    r = 0;

  float delta_eta_mm = r * cos(alpha);
  float delta_phi_mm = r * sin(alpha);

  // Particles with negative eta are expected to have the same shape as those
  // with positive eta after transformation: delta_eta --> -delta_eta
  if (center_eta < 0.)
    delta_eta_mm = -delta_eta_mm;
  // We derive the shower shapes for electrons and positively charged hadrons.
  // Particle with the opposite charge are expected to have the same shower
  // shape after the transformation: delta_phi --> -delta_phi
  if ((charge < 0. && pdgId != 11) || pdgId == -11)
    delta_phi_mm = -delta_phi_mm;

  const float dist000 = TMath::Sqrt(center_r * center_r + center_z * center_z);
  const float eta_jakobi = TMath::Abs(2.0 * TMath::Exp(-center_eta)
                                      / (1.0 + TMath::Exp(-2 * center_eta)));

  const float delta_eta = delta_eta_mm / eta_jakobi / dist000;
  const float delta_phi = delta_phi_mm / center_r;

  hit.setEtaPhiZE(
      center_eta + delta_eta, center_phi + delta_phi, center_z, hit.E());

  FCS_MSG_DEBUG("HIT: E=" << hit.E() << " cs=" << cs << " eta=" << hit.eta()
                          << " phi=" << hit.phi() << " z=" << hit.z()
                          << " r=" << r << " alpha=" << alpha);

  return FCSSuccess;
}

bool TFCSHistoLateralShapeParametrization::Initialize(TH2* hist)
{
  if (!hist)
    return false;
  m_hist.Initialize(hist);
  if (m_hist.get_HistoContents().empty())
    return false;

  set_number_of_hits(hist->Integral());

  return true;
}

bool TFCSHistoLateralShapeParametrization::Initialize(const char* filepath,
                                                      const char* histname)
{
  // input file with histogram to fit
  std::unique_ptr<TFile> inputfile(TFile::Open(filepath, "READ"));
  if (inputfile == nullptr)
    return false;

  // histogram with hit pattern
  TH2* inputShape = (TH2*)inputfile->Get(histname);
  if (inputShape == nullptr)
    return false;

  bool OK = Initialize(inputShape);

  inputfile->Close();

  return OK;
}

void TFCSHistoLateralShapeParametrization::Print(Option_t* option) const
{
  TString opt(option);
  bool shortprint = opt.Index("short") >= 0;
  bool longprint =
      msgLvl(FCS_MSG::DEBUG) || (msgLvl(FCS_MSG::INFO) && !shortprint);
  TString optprint = opt;
  optprint.ReplaceAll("short", "");
  TFCSLateralShapeParametrizationHitBase::Print(option);

  if (longprint) {
    if (is_phi_symmetric()) {
      FCS_MSG_INFO(optprint
                   << "  Histo: " << m_hist.get_HistoBordersx().size() - 1
                   << "*" << m_hist.get_HistoBordersy().size() - 1
                   << " bins, #hits=" << m_nhits << ", r scale=" << m_r_scale
                   << ", r offset=" << m_r_offset << "mm (phi symmetric)");
    } else {
      FCS_MSG_INFO(optprint
                   << "  Histo: " << m_hist.get_HistoBordersx().size() - 1
                   << "*" << m_hist.get_HistoBordersy().size() - 1
                   << " bins, #hits=" << m_nhits << ", r scale=" << m_r_scale
                   << ", r offset=" << m_r_offset << "mm (not phi symmetric)");
    }
  }
}
