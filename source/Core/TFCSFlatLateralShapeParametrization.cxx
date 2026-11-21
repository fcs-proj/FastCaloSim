// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#include <cmath>
#include <numbers>

#include "FastCaloSim/Core/TFCSFlatLateralShapeParametrization.h"

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandPoisson.h"
#include "FastCaloSim/Core/TFCSExtrapolationState.h"
#include "FastCaloSim/Core/TFCSSimulationState.h"
#include "TFile.h"
#include "TH2.h"

//=============================================
//======= TFCSFlatLateralShapeParametrization =========
//=============================================

TFCSFlatLateralShapeParametrization::TFCSFlatLateralShapeParametrization(
    const char* name, const char* title)
    : TFCSLateralShapeParametrizationHitBase(name, title)
    , m_nhits(0)
    , m_scale(1)
{
}

TFCSFlatLateralShapeParametrization::~TFCSFlatLateralShapeParametrization() =
    default;

auto TFCSFlatLateralShapeParametrization::get_number_of_hits(
    TFCSSimulationState& simulstate,
    const TFCSTruthState* /*truth*/,
    const TFCSExtrapolationState* /*extrapol*/) const -> int
{
  if (!simulstate.randomEngine()) {
    return -1;
  }

  return CLHEP::RandPoisson::shoot(simulstate.randomEngine(), m_nhits);
}

void TFCSFlatLateralShapeParametrization::set_number_of_hits(float nhits)
{
  m_nhits = nhits;
}

void TFCSFlatLateralShapeParametrization::set_dR(float _dR)
{
  m_dR = _dR;
}

void TFCSFlatLateralShapeParametrization::set_scale(float _scale)
{
  m_scale = _scale;
}

auto TFCSFlatLateralShapeParametrization::simulate_hit(
    Hit& hit,
    TFCSSimulationState& simulstate,
    const TFCSTruthState* /*truth*/,
    const TFCSExtrapolationState* /*extrapol*/) -> FCSReturnCode
{
  if (!simulstate.randomEngine()) {
    return FCSFatal;
  }

  const int cs = calosample();
  const double center_eta = hit.center_eta();
  const double center_phi = hit.center_phi();
  const double center_r = hit.center_r();
  const double center_z = hit.center_z();

  if (std::isnan(center_r) or std::isnan(center_z) or std::isnan(center_eta)
      or std::isnan(center_phi))
  {  // Check if extrapolation fails
    return FCSFatal;
  }

  float alpha, r;

  alpha =
      2 * std::numbers::pi * CLHEP::RandFlat::shoot(simulstate.randomEngine());
  r = m_dR * CLHEP::RandFlat::shoot(simulstate.randomEngine());

  float delta_eta = r * std::cos(alpha);
  float delta_phi = r * std::sin(alpha);

  hit.setEtaPhiZE(center_eta + delta_eta,
                  center_phi + delta_phi,
                  center_z,
                  hit.E() * m_scale);

  FCS_MSG_DEBUG("HIT: E=" << hit.E() << " cs=" << cs << " eta=" << hit.eta()
                          << " phi=" << hit.phi() << " z=" << hit.z()
                          << " r=" << r << " alpha=" << alpha);

  return FCSSuccess;
}

void TFCSFlatLateralShapeParametrization::Print(Option_t* option) const
{
  const auto [shortprint, longprint, optprint] = print_options(option);
  TFCSLateralShapeParametrizationHitBase::Print(option);

  if (longprint) {
    FCS_MSG_INFO(optprint << "  dR=" << m_dR << " scale factor=" << m_scale
                          << ", #hits=" << m_nhits);
  }
}
