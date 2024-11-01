// Copyright (c) 2024 CERN for the benefit of the FastCaloSim project

#include "FastCaloSim/Core/TFCS2DFunctionLateralShapeParametrization.h"

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandPoisson.h"
#include "FastCaloSim/Core/TFCSExtrapolationState.h"
#include "FastCaloSim/Core/TFCSSimulationState.h"
#include "FastCaloSim/Definitions/ParticleData.h"
#include "TFile.h"
#include "TH2.h"
#include "TMath.h"

//=============================================
//======= TFCS2DFunctionLateralShapeParametrization =========
//=============================================

TFCS2DFunctionLateralShapeParametrization::
    TFCS2DFunctionLateralShapeParametrization(const char* name,
                                              const char* title)
    : TFCSLateralShapeParametrizationHitBase(name, title)
    , m_function(nullptr)
    , m_nhits(0)
{
  TFCS2DFunctionLateralShapeParametrization::reset_phi_symmetric();
}

TFCS2DFunctionLateralShapeParametrization::
    ~TFCS2DFunctionLateralShapeParametrization()
{
  if (m_function)
    delete m_function;
  m_function = nullptr;
}

double TFCS2DFunctionLateralShapeParametrization::get_sigma2_fluctuation(
    TFCSSimulationState& /*simulstate*/,
    const TFCSTruthState* /*truth*/,
    const TFCSExtrapolationState* /*extrapol*/) const
{
  return 1.0 / m_nhits;
}

int TFCS2DFunctionLateralShapeParametrization::get_number_of_hits(
    TFCSSimulationState& simulstate,
    const TFCSTruthState* /*truth*/,
    const TFCSExtrapolationState* /*extrapol*/) const
{
  if (!simulstate.randomEngine()) {
    return -1;
  }

  return CLHEP::RandPoisson::shoot(simulstate.randomEngine(), m_nhits);
}

void TFCS2DFunctionLateralShapeParametrization::set_number_of_hits(float nhits)
{
  m_nhits = nhits;
}

FCSReturnCode TFCS2DFunctionLateralShapeParametrization::simulate_hit(
    Hit& hit,
    TFCSSimulationState& simulstate,
    const TFCSTruthState* truth,
    const TFCSExtrapolationState* /*extrapol*/)
{
  if (!simulstate.randomEngine()) {
    return FCSFatal;
  }
  if (!m_function) {
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
      m_function->rnd_to_fct(alpha, r, rnd1, rnd2);
      alpha = -alpha;
    } else {  // Fill positive phi half of shape
      rnd2 *= 2;
      m_function->rnd_to_fct(alpha, r, rnd1, rnd2);
    }
  } else {
    m_function->rnd_to_fct(alpha, r, rnd1, rnd2);
  }
  if (TMath::IsNaN(alpha) || TMath::IsNaN(r)) {
    FCS_MSG_ERROR("  2D function, #hits=" << m_nhits << " alpha=" << alpha
                                          << " r=" << r << " rnd1=" << rnd1
                                          << " rnd2=" << rnd2);
    alpha = 0;
    r = 0.001;

    FCS_MSG_ERROR("  This error could probably be retried");
    return FCSFatal;
  }

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

bool TFCS2DFunctionLateralShapeParametrization::Initialize(TFCS2DFunction* func,
                                                           float nhits)
{
  if (!func)
    return false;
  if (m_function)
    delete m_function;
  m_function = func;

  if (nhits > 0)
    set_number_of_hits(nhits);

  return true;
}

void TFCS2DFunctionLateralShapeParametrization::Print(Option_t* option) const
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
      FCS_MSG_INFO(optprint << "  2D function, #hits=" << m_nhits
                            << " (phi symmetric)");
    } else {
      FCS_MSG_INFO(optprint << "  2D function, #hits=" << m_nhits
                            << " (not phi symmetric)");
    }
  }
}
