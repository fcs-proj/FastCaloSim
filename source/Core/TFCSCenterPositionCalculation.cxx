/*
  Copyright (C) 2002-2018 CERN for the benefit of the ATLAS collaboration
*/

#include "FastCaloSim/Core/TFCSCenterPositionCalculation.h"

#include "FastCaloSim/Core/TFCSExtrapolationState.h"
#include "FastCaloSim/Core/TFCSSimulationState.h"
#include "FastCaloSim/Geometry/Cell.h"

//=============================================
//======= TFCSCenterPositionCalculation =========
//=============================================

TFCSCenterPositionCalculation::TFCSCenterPositionCalculation(const char* name,
                                                             const char* title)
    : TFCSLateralShapeParametrizationHitBase(name, title)
    , m_extrapWeight(0.5)
{
}

FCSReturnCode TFCSCenterPositionCalculation::simulate_hit(
    Hit& hit,
    TFCSSimulationState& /*simulstate*/,
    const TFCSTruthState* /*truth*/,
    const TFCSExtrapolationState* extrapol)
{
  const int cs = calosample();

  double r = (1. - m_extrapWeight) * extrapol->r(cs, Cell::SubPos::ENT)
      + m_extrapWeight * extrapol->r(cs, Cell::SubPos::EXT);
  double z = (1. - m_extrapWeight) * extrapol->z(cs, Cell::SubPos::ENT)
      + m_extrapWeight * extrapol->z(cs, Cell::SubPos::EXT);
  double eta = (1. - m_extrapWeight) * extrapol->eta(cs, Cell::SubPos::ENT)
      + m_extrapWeight * extrapol->eta(cs, Cell::SubPos::EXT);
  double phi = (1. - m_extrapWeight) * extrapol->phi(cs, Cell::SubPos::ENT)
      + m_extrapWeight * extrapol->phi(cs, Cell::SubPos::EXT);

  if (!std::isfinite(r) || !std::isfinite(z) || !std::isfinite(eta)
      || !std::isfinite(phi))
  {
    ATH_MSG_WARNING(
        "Extrapolator contains NaN or infinite number.\nSetting "
        "center position to calo boundary.");
    ATH_MSG_WARNING("Before fix: center_r: "
                    << r << " center_z: " << z << " center_phi: " << phi
                    << " center_eta: " << eta << " weight: " << m_extrapWeight
                    << " cs: " << cs);
    // If extrapolator fails we can set position to calo boundary
    r = extrapol->IDCaloBoundary_r();
    z = extrapol->IDCaloBoundary_z();
    eta = extrapol->IDCaloBoundary_eta();
    phi = extrapol->IDCaloBoundary_phi();

    ATH_MSG_WARNING("After fix: center_r: "
                    << r << " center_z: " << z << " center_phi: " << phi
                    << " center_eta: " << eta << " weight: " << m_extrapWeight
                    << " cs: " << cs);
  }

  hit.setCenter_r(r);
  hit.setCenter_z(z);
  hit.setCenter_eta(eta);
  hit.setCenter_phi(phi);

  ATH_MSG_DEBUG("TFCSCenterPositionCalculation: center_r: "
                << hit.center_r() << " center_z: " << hit.center_z()
                << " center_phi: " << hit.center_phi()
                << " center_eta: " << hit.center_eta()
                << " weight: " << m_extrapWeight << " cs: " << cs);

  return FCSSuccess;
}

void TFCSCenterPositionCalculation::Print(Option_t* option) const
{
  TString opt(option);
  bool shortprint = opt.Index("short") >= 0;
  bool longprint = msgLvl(MSG::DEBUG) || (msgLvl(MSG::INFO) && !shortprint);
  TString optprint = opt;
  optprint.ReplaceAll("short", "");
  TFCSLateralShapeParametrizationHitBase::Print(option);

  if (longprint)
    ATH_MSG_INFO(optprint << "  Weight for extrapolated position: "
                          << m_extrapWeight);
}
