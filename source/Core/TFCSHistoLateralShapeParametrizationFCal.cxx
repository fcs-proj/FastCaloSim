// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#include <cmath>

#include "FastCaloSim/Core/TFCSHistoLateralShapeParametrizationFCal.h"

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandPoisson.h"
#include "FastCaloSim/Core/TFCSExtrapolationState.h"
#include "FastCaloSim/Core/TFCSSimulationState.h"
#include "FastCaloSim/Definitions/ParticleData.h"

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
    ~TFCSHistoLateralShapeParametrizationFCal() = default;

auto TFCSHistoLateralShapeParametrizationFCal::simulate_hit(
    Hit& hit,
    TFCSSimulationState& simulstate,
    const TFCSTruthState* truth,
    const TFCSExtrapolationState* /*extrapol*/) -> FCSReturnCode
{
  if (simulstate.randomEngine() == nullptr) {
    return FCSFatal;
  }

  const int pdgId = truth->pdgid();
  const double charge = ParticleData::charge(pdgId);

  const std::size_t cs = calosample();
  const double center_phi = hit.center_phi();
  const double center_r = hit.center_r();
  const double center_z = hit.center_z();

  // Set by m_hist.rnd_to_fct()
  double alpha = 0;
  double r = 0;
  const double rnd0 = CLHEP::RandFlat::shoot(simulstate.randomEngine());
  double rnd1 = CLHEP::RandFlat::shoot(simulstate.randomEngine());
  if (is_phi_symmetric()) {
    if (rnd1 >= 0.5) {  // Fill negative phi half of shape
      rnd1 -= 0.5;
      rnd1 *= 2;
      m_hist.rnd_to_fct(alpha, r, rnd0, rnd1);
      alpha = -alpha;
    } else {  // Fill positive phi half of shape
      rnd1 *= 2;
      m_hist.rnd_to_fct(alpha, r, rnd0, rnd1);
    }
  } else {
    m_hist.rnd_to_fct(alpha, r, rnd0, rnd1);
  }
  if (std::isnan(alpha) || std::isnan(r)) {
    FCS_MSG_ERROR("  Histogram: "
                  << m_hist.get_HistoBordersx().size() - 1 << "*"
                  << m_hist.get_HistoBordersy().size() - 1
                  << " bins, #hits=" << m_nhits << " alpha=" << alpha
                  << " r=" << r << " rnd0=" << rnd0 << " rnd1=" << rnd1);

    FCS_MSG_ERROR("  This error could probably be retried");
    return FCSFatal;
  }

  const double hit_r = r * std::cos(alpha) + center_r;
  double delta_phi = r * std::sin(alpha) / center_r;
  // We derive the shower shapes for electrons and positively charged hadrons.
  // Particle with the opposite charge are expected to have the same shower
  // shape after the transformation: delta_phi --> -delta_phi
  if ((charge < 0. && pdgId != 11) || pdgId == -11) {
    delta_phi = -delta_phi;
  }
  const double hit_phi = delta_phi + center_phi;

  hit.setXYZE(
      hit_r * std::cos(hit_phi), hit_r * std::sin(hit_phi), center_z, hit.E());

  FCS_MSG_DEBUG("HIT: E=" << hit.E() << " cs=" << cs << " x=" << hit.x()
                          << " y=" << hit.y() << " z=" << hit.z() << " r=" << r
                          << " alpha=" << alpha);

  return FCSSuccess;
}
