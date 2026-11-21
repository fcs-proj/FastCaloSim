// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#include <cmath>

#include "FastCaloSim/Core/TFCS2DFunctionLateralShapeParametrization.h"

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandPoisson.h"
#include "FastCaloSim/Core/TFCSExtrapolationState.h"
#include "FastCaloSim/Core/TFCSSimulationState.h"
#include "FastCaloSim/Definitions/ParticleData.h"
#include "TFile.h"
#include "TH2.h"

//=============================================
//======= TFCS2DFunctionLateralShapeParametrization =========
//=============================================

TFCS2DFunctionLateralShapeParametrization::
    TFCS2DFunctionLateralShapeParametrization(const char* name,
                                              const char* title)
    : TFCSLateralShapeParametrizationHitBase(name, title)

{
  TFCS2DFunctionLateralShapeParametrization::reset_phi_symmetric();
}

TFCS2DFunctionLateralShapeParametrization::
    ~TFCS2DFunctionLateralShapeParametrization()
{
  {
    delete m_function;
  }
  m_function = nullptr;
}

auto TFCS2DFunctionLateralShapeParametrization::get_sigma2_fluctuation(
    TFCSSimulationState& /*simulstate*/,
    const TFCSTruthState* /*truth*/,
    const TFCSExtrapolationState* /*extrapol*/) const -> double
{
  return 1.0 / static_cast<double>(m_nhits);
}

auto TFCS2DFunctionLateralShapeParametrization::get_number_of_hits(
    TFCSSimulationState& simulstate,
    const TFCSTruthState* /*truth*/,
    const TFCSExtrapolationState* /*extrapol*/) const -> int
{
  if (simulstate.randomEngine() == nullptr) {
    return -1;
  }

  return static_cast<int>(CLHEP::RandPoisson::shoot(
      simulstate.randomEngine(), static_cast<double>(m_nhits)));
}

void TFCS2DFunctionLateralShapeParametrization::set_number_of_hits(float nhits)
{
  m_nhits = nhits;
}

auto TFCS2DFunctionLateralShapeParametrization::simulate_hit(
    Hit& hit,
    TFCSSimulationState& simulstate,
    const TFCSTruthState* truth,
    const TFCSExtrapolationState* /*extrapol*/) -> FCSReturnCode
{
  if (simulstate.randomEngine() == nullptr) {
    return FCSFatal;
  }
  if (m_function == nullptr) {
    return FCSFatal;
  }

  const int pdgId = truth->pdgid();
  const double charge = ParticleData::charge(pdgId);

  const std::size_t cs = calosample();
  const double center_eta = hit.center_eta();
  const double center_phi = hit.center_phi();
  const double center_r = hit.center_r();
  const double center_z = hit.center_z();

  if (std::isnan(center_r) or std::isnan(center_z) or std::isnan(center_eta)
      or std::isnan(center_phi))
  {  // Check if extrapolation fails
    return FCSFatal;
  }

  // Set by m_function->rnd_to_fct()
  double alpha = 0;
  double r = 0;
  const double rnd0 = CLHEP::RandFlat::shoot(simulstate.randomEngine());
  double rnd1 = CLHEP::RandFlat::shoot(simulstate.randomEngine());
  if (is_phi_symmetric()) {
    if (rnd1 >= 0.5) {  // Fill negative phi half of shape
      rnd1 -= 0.5;
      rnd1 *= 2;
      m_function->rnd_to_fct(alpha, r, rnd0, rnd1);
      alpha = -alpha;
    } else {  // Fill positive phi half of shape
      rnd1 *= 2;
      m_function->rnd_to_fct(alpha, r, rnd0, rnd1);
    }
  } else {
    m_function->rnd_to_fct(alpha, r, rnd0, rnd1);
  }
  if (std::isnan(alpha) || std::isnan(r)) {
    FCS_MSG_ERROR("  2D function, #hits=" << m_nhits << " alpha=" << alpha
                                          << " r=" << r << " rnd0=" << rnd0
                                          << " rnd1=" << rnd1);

    FCS_MSG_ERROR("  This error could probably be retried");
    return FCSFatal;
  }

  double delta_eta_mm = r * std::cos(alpha);
  double delta_phi_mm = r * std::sin(alpha);

  // Particles with negative eta are expected to have the same shape as those
  // with positive eta after transformation: delta_eta --> -delta_eta
  if (center_eta < 0.) {
    delta_eta_mm = -delta_eta_mm;
  }
  // We derive the shower shapes for electrons and positively charged hadrons.
  // Particle with the opposite charge are expected to have the same shower
  // shape after the transformation: delta_phi --> -delta_phi
  if ((charge < 0. && pdgId != 11) || pdgId == -11) {
    delta_phi_mm = -delta_phi_mm;
  }

  const double dist000 = FastCaloSim::dist000(center_r, center_z);
  const double eta_jakobi = FastCaloSim::eta_jakobi(center_eta);

  const double delta_eta = delta_eta_mm / eta_jakobi / dist000;
  const double delta_phi = delta_phi_mm / center_r;

  hit.setEtaPhiZE(
      center_eta + delta_eta, center_phi + delta_phi, center_z, hit.E());

  FCS_MSG_DEBUG("HIT: E=" << hit.E() << " cs=" << cs << " eta=" << hit.eta()
                          << " phi=" << hit.phi() << " z=" << hit.z()
                          << " r=" << r << " alpha=" << alpha);

  return FCSSuccess;
}

auto TFCS2DFunctionLateralShapeParametrization::Initialize(TFCS2DFunction* func,
                                                           float nhits) -> bool
{
  if (func == nullptr) {
    return false;
  }
  {
    delete m_function;
  }
  m_function = func;

  if (nhits > 0) {
    set_number_of_hits(nhits);
  }

  return true;
}

void TFCS2DFunctionLateralShapeParametrization::Print(Option_t* option) const
{
  const auto [shortprint, longprint, optprint] = print_options(option);
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
