// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

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

auto TFCSHitCellMapping::simulate_hit(
    Hit& hit,
    TFCSSimulationState& simulstate,
    const TFCSTruthState* /*truth*/,
    const TFCSExtrapolationState* /*extrapol*/) -> FCSReturnCode
{
  FCS_MSG_DEBUG("Got hit with E=" << hit.E() << " eta=" << hit.eta()
                                  << " phi=" << hit.phi());

  // Position where we will perform the lookup
  Position const lookup_pos {.m_x = 0,
                             .m_y = 0,
                             .m_z = 0,
                             .m_eta = hit.eta(),
                             .m_phi = hit.phi(),
                             .m_r = 0};

  // Get the best matching cell
  const auto& cell = m_geo->get_cell(calosample(), lookup_pos);
  FCS_MSG_DEBUG(cell);

  // Get hit-cell boundary proximity
  // < 0 means we are inside the cell
  // > 0 means we are outside the cell
  double const proximity = cell.boundary_proximity(lookup_pos);

  FCS_MSG_DEBUG("Hit-cell distance in eta-phi is: " << proximity);

  // If the distance is positive then we are using the nearest cell rather
  // than are inside a cell If we are more than 0.005mm from the nearest cell
  // we don't create a hit to avoid the build-up of energy in edge cells
  // For FCSV2 another hit can be created but with a cutoff to avoid looping,
  // for FastCaloGAN the rest of the hits in the layer will be scaled up by the
  // energy renormalization step.
  if (proximity < 0.005) {
    // Cell energies are stored in single precision
    simulstate.deposit(cell.id(), hit.E());
  } else {
    hit.setXYZE(hit.x(), hit.y(), hit.z(), 0.0);
  }
  return FCSSuccess;
}

auto TFCSHitCellMapping::operator==(const TFCSParametrizationBase& ref) const
    -> bool
{
  if (TFCSParametrizationBase::compare(ref)) {
    return true;
  }
  if (!TFCSParametrization::compare(ref)) {
    return false;
  }
  return TFCSLateralShapeParametrization::compare(ref);
}

void TFCSHitCellMapping::Print(Option_t* option) const
{
  const auto [shortprint, longprint, optprint] = print_options(option);
  TFCSLateralShapeParametrizationHitBase::Print(option);

  if (longprint) {
    FCS_MSG_INFO(optprint << "  geo=" << m_geo);
  }
}
