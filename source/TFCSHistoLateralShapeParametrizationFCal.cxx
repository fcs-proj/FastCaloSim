/*
  Copyright (C) 2002-2018 CERN for the benefit of the ATLAS collaboration
*/

#include "FastCaloSim/TFCSHistoLateralShapeParametrizationFCal.h"

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandPoisson.h"
#include "FastCaloSim/TFCSExtrapolationState.h"
#include "FastCaloSim/TFCSSimulationState.h"
#include "HepPDT/ParticleData.hh"
#include "HepPDT/ParticleDataTable.hh"
#include "TMath.h"

//=============================================
//======= TFCSHistoLateralShapeParametrizationFCal =========
//=============================================

TFCSHistoLateralShapeParametrizationFCal::
    TFCSHistoLateralShapeParametrizationFCal(const char* name,
                                             const char* title)
    : TFCSHistoLateralShapeParametrization(name, title)
{
}

TFCSHistoLateralShapeParametrizationFCal::
    ~TFCSHistoLateralShapeParametrizationFCal()
{
}

FCSReturnCode TFCSHistoLateralShapeParametrizationFCal::simulate_hit(
    Hit& hit,
    TFCSSimulationState& simulstate,
    const TFCSTruthState* truth,
    const TFCSExtrapolationState* /*extrapol*/)
{
  if (!simulstate.randomEngine()) {
    return FCSFatal;
  }

  const int pdgId = truth->pdgid();
  const double charge = HepPDT::ParticleID(pdgId).charge();

  const int cs = calosample();
  // const double center_phi=0.5*( extrapol->phi(cs, CaloSubPos::SUBPOS_ENT) +
  // extrapol->phi(cs, CaloSubPos::SUBPOS_EXT) ); const double center_r=0.5*(
  // extrapol->r(cs, CaloSubPos::SUBPOS_ENT) + extrapol->r(cs,
  // CaloSubPos::SUBPOS_EXT) ); const double center_z=0.5*( extrapol->z(cs,
  // CaloSubPos::SUBPOS_ENT) + extrapol->z(cs, CaloSubPos::SUBPOS_EXT) );
  const double center_phi = hit.center_phi();
  const double center_r = hit.center_r();
  const double center_z = hit.center_z();

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
    ATH_MSG_ERROR("  Histogram: "
                  << m_hist.get_HistoBordersx().size() - 1 << "*"
                  << m_hist.get_HistoBordersy().size() - 1
                  << " bins, #hits=" << m_nhits << " alpha=" << alpha
                  << " r=" << r << " rnd1=" << rnd1 << " rnd2=" << rnd2);
    alpha = 0;
    r = 0.001;

    ATH_MSG_ERROR("  This error could probably be retried");
    return FCSFatal;
  }

  const float hit_r = r * cos(alpha) + center_r;
  float delta_phi = r * sin(alpha) / center_r;
  // We derive the shower shapes for electrons and positively charged hadrons.
  // Particle with the opposite charge are expected to have the same shower
  // shape after the transformation: delta_phi --> -delta_phi
  if ((charge < 0. && pdgId != 11) || pdgId == -11)
    delta_phi = -delta_phi;
  const float hit_phi = delta_phi + center_phi;

  hit.setXYZE(hit_r * cos(hit_phi), hit_r * sin(hit_phi), center_z, hit.E());

  ATH_MSG_DEBUG("HIT: E=" << hit.E() << " cs=" << cs << " x=" << hit.x()
                          << " y=" << hit.y() << " z=" << hit.z() << " r=" << r
                          << " alpha=" << alpha);

  return FCSSuccess;
}
