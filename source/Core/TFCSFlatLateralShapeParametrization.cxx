/*
  Copyright (C) 2002-2019 CERN for the benefit of the ATLAS collaboration
*/

#include "FastCaloSim/Core/TFCSFlatLateralShapeParametrization.h"

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandPoisson.h"
#include "FastCaloSim/Geometry/FastCaloSim_CaloCell_ID.h"
#include "FastCaloSim/Core/TFCSExtrapolationState.h"
#include "FastCaloSim/Core/TFCSSimulationState.h"
#include "TFile.h"
#include "TH2.h"
#include "TMath.h"

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

TFCSFlatLateralShapeParametrization::~TFCSFlatLateralShapeParametrization() {}

int TFCSFlatLateralShapeParametrization::get_number_of_hits(
    TFCSSimulationState& simulstate,
    const TFCSTruthState* /*truth*/,
    const TFCSExtrapolationState* /*extrapol*/) const
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

FCSReturnCode TFCSFlatLateralShapeParametrization::simulate_hit(
    Hit& hit,
    TFCSSimulationState& simulstate,
    const TFCSTruthState* /*truth*/,
    const TFCSExtrapolationState* /*extrapol*/)
{
  if (!simulstate.randomEngine()) {
    return FCSFatal;
  }

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

  float alpha, r;

  alpha = 2 * TMath::Pi() * CLHEP::RandFlat::shoot(simulstate.randomEngine());
  r = m_dR * CLHEP::RandFlat::shoot(simulstate.randomEngine());

  float delta_eta = r * cos(alpha);
  float delta_phi = r * sin(alpha);

  hit.setEtaPhiZE(center_eta + delta_eta,
                  center_phi + delta_phi,
                  center_z,
                  hit.E() * m_scale);

  ATH_MSG_DEBUG("HIT: E=" << hit.E() << " cs=" << cs << " eta=" << hit.eta()
                          << " phi=" << hit.phi() << " z=" << hit.z()
                          << " r=" << r << " alpha=" << alpha);

  return FCSSuccess;
}

void TFCSFlatLateralShapeParametrization::Print(Option_t* option) const
{
  TString opt(option);
  bool shortprint = opt.Index("short") >= 0;
  bool longprint = msgLvl(MSG::DEBUG) || (msgLvl(MSG::INFO) && !shortprint);
  TString optprint = opt;
  optprint.ReplaceAll("short", "");
  TFCSLateralShapeParametrizationHitBase::Print(option);

  if (longprint) {
    ATH_MSG_INFO(optprint << "  dR=" << m_dR << " scale factor=" << m_scale
                          << ", #hits=" << m_nhits);
  }
}
