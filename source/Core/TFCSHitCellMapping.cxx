#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-equal"

/*
  Copyright (C) 2002-2023 CERN for the benefit of the ATLAS collaboration
*/

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
  ATH_MSG_DEBUG("Got hit with E=" << hit.E() << " eta=" << hit.eta()
                                  << " phi=" << hit.phi());

  // Get the best matching cell
  const auto& cell = m_geo->get_cell(calosample(), hit);
  ATH_MSG_DEBUG(cell);

  // Get hit-cell distance
  // < 0 means we are inside the cell
  // > 0 means we are outside the cell
  double distance = dist(cell, hit);

  ATH_MSG_DEBUG("Hit-cell distance in eta-phi is: " << distance << " mm");

  // If the distance is positive then we are using the nearest cell rather
  // than are inside a cell If we are more than 0.005mm from the nearest cell
  // we don't create a hit to avoid the build-up of energy in edge cells
  // For FCSV2 another hit can be created but with a cutoff to avoid looping,
  // for FastCaloGAN the rest of the hits in the layer will be scaled up by the
  // energy renormalization step.
  if (distance < 0.005) {
    simulstate.deposit(cell.id(), hit.E());
  } else {
    hit.setXYZE(hit.x(), hit.y(), hit.z(), 0.0);
  }
  return FCSSuccess;
}

auto TFCSHitCellMapping::mod_diff(double angle1, double angle2) -> double
{
  double difference = std::fmod(std::abs(angle1 - angle2), 2 * M_PI);
  return difference > M_PI ? 2 * M_PI - difference : difference;
}

auto TFCSHitCellMapping::dist(const Cell& cell, const Hit& hit) -> double
{
  // Calculate the distances in eta and phi directions
  double delta_eta = std::abs(hit.eta() - cell.eta());
  double delta_phi = mod_diff(hit.phi(), cell.phi());

  // Closest distance to a boundary
  return std::min(delta_eta - cell.deta(), delta_phi - cell.dphi());
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
