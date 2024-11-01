#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-equal"

// Copyright (c) 2024 CERN for the benefit of the FastCaloSim project

#include "FastCaloSim/Core/TFCSHitCellMapping.h"

#include "FastCaloSim/Core/TFCSSimulationState.h"
#include "FastCaloSim/Geometry/CaloGeo.h"

//=============================================
//======= TFCSHitCellMapping =========
//=============================================

TFCSHitCellMapping::TFCSHitCellMapping(const char* name,
                                       const char* title,
                                       CaloGeo* geo)
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
  FCS_MSG_DEBUG("Got hit with E=" << hit.E() << " eta=" << hit.eta()
                                  << " phi=" << hit.phi());

  // Position where we will perform the lookup
  Position lookup_pos {0, 0, 0, hit.eta(), hit.phi(), 0};

  // Get the best matching cell
  const auto& cell = m_geo->get_cell(calosample(), lookup_pos);
  FCS_MSG_DEBUG(cell);

  // Get hit-cell boundary proximity
  // < 0 means we are inside the cell
  // > 0 means we are outside the cell
  double proximity = cell.boundary_proximity(lookup_pos);

  FCS_MSG_DEBUG("Hit-cell distance in eta-phi is: " << proximity);

  // If the distance is positive then we are using the nearest cell rather
  // than are inside a cell If we are more than 0.005mm from the nearest cell
  // we don't create a hit to avoid the build-up of energy in edge cells
  // For FCSV2 another hit can be created but with a cutoff to avoid looping,
  // for FastCaloGAN the rest of the hits in the layer will be scaled up by the
  // energy renormalization step.
  if (proximity < 0.005) {
    simulstate.deposit(cell.id(), hit.E());
  } else {
    hit.setXYZE(hit.x(), hit.y(), hit.z(), 0.0);
  }
  return FCSSuccess;
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
  bool longprint =
      msgLvl(FCS_MSG::DEBUG) || (msgLvl(FCS_MSG::INFO) && !shortprint);
  TString optprint = opt;
  optprint.ReplaceAll("short", "");
  TFCSLateralShapeParametrizationHitBase::Print(option);

  if (longprint)
    FCS_MSG_INFO(optprint << "  geo=" << m_geo);
}

#pragma GCC diagnostic push
