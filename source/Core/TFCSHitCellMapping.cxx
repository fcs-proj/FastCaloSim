#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-equal"

/*
  Copyright (C) 2002-2023 CERN for the benefit of the ATLAS collaboration
*/

#include "FastCaloSim/Core/TFCSHitCellMapping.h"

#include "FastCaloSim/Core/TFCSSimulationState.h"
#include "FastCaloSim/Geometry/ICaloGeometry.h"

//=============================================
//======= TFCSHitCellMapping =========
//=============================================

TFCSHitCellMapping::TFCSHitCellMapping(const char* name,
                                       const char* title,
                                       ICaloGeometry* geo)
    : TFCSLateralShapeParametrizationHitBase(name, title)
    , m_geo(geo)
{
  set_match_all_pdgid();
}

FCSReturnCode TFCSHitCellMapping::simulate_hit(
    Hit& hit,
    TFCSSimulationState& simulstate,
    const TFCSTruthState* /*truth*/,
    const TFCSExtrapolationState* /*extrapol*/)
{
  int cs = calosample();
  float distance;
  const CaloDetDescrElement* cellele =
      m_geo->getDDE(cs, hit.eta(), hit.phi(), &distance);
  ATH_MSG_DEBUG("HIT: cellele=" << cellele << " E=" << hit.E() << " cs=" << cs
                                << " eta=" << hit.eta()
                                << " phi=" << hit.phi());
  if (cellele) {
    // If the distance is positive then we are using the nearest cell rather
    // than are inside a cell If we are more than 0.005mm from the nearest cell
    // we don't create a hit to avoid the build-up of energy in edge cells For
    // FCSV2 another hit can be created but with a cutoff to avoid looping, for
    // FastCaloGAN the rest of the hits in the layer will be scaled up by the
    // energy renormalization step.
    if (distance < 0.005) {
      simulstate.deposit(cellele, hit.E());
    } else {
      hit.setXYZE(hit.x(), hit.y(), hit.z(), 0.0);
    }
    return FCSSuccess;
  } else {
    ATH_MSG_ERROR(
        "TFCSLateralShapeParametrizationHitCellMapping::simulate_hit: cellele="
        << cellele << " E=" << hit.E() << " cs=" << cs << " eta=" << hit.eta()
        << " phi=" << hit.phi());
    return FCSFatal;
  }
}

bool TFCSHitCellMapping::operator==(const TFCSParametrizationBase& ref) const
{
  if (TFCSParametrizationBase::compare(ref))
    return true;
  if (!TFCSParametrization::compare(ref))
    return false;
  if (!TFCSLateralShapeParametrization::compare(ref))
    return false;

  return true;
}

void TFCSHitCellMapping::Print(Option_t* option) const
{
  TString opt(option);
  bool shortprint = opt.Index("short") >= 0;
  bool longprint = msgLvl(MSG::DEBUG) || (msgLvl(MSG::INFO) && !shortprint);
  TString optprint = opt;
  optprint.ReplaceAll("short", "");
  TFCSLateralShapeParametrizationHitBase::Print(option);

  if (longprint)
    ATH_MSG_INFO(optprint << "  geo=" << m_geo);
}

#pragma GCC diagnostic push
