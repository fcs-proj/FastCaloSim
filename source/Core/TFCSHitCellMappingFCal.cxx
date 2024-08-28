#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-equal"

/*
  Copyright (C) 2002-2018 CERN for the benefit of the ATLAS collaboration
*/

#include "FastCaloSim/Core/TFCSHitCellMappingFCal.h"

#include "FastCaloSim/Core/TFCSSimulationState.h"
#include "FastCaloSim/Geometry/CaloGeo.h"

//=============================================
//======= TFCSHitCellMappingFCal =========
//=============================================

FCSReturnCode TFCSHitCellMappingFCal::simulate_hit(
    Hit& hit,
    TFCSSimulationState& simulstate,
    const TFCSTruthState* /*truth*/,
    const TFCSExtrapolationState* /*extrapol*/)
{
  ATH_MSG_DEBUG("Got hit with E=" << hit.E() << " x=" << hit.x()
                                  << " y=" << hit.y());

  // Get the best matching cell
  // TODO: implement proper cell matching for FCAL
  const auto& cell = m_geo->get_cell(calosample(), hit);
  ATH_MSG_DEBUG(cell);

  // Get hit-cell boundary proximity
  // < 0 means we are inside the cell
  // > 0 means we are outside the cell
  double proximity = cell.boundary_proximity(hit);

  ATH_MSG_DEBUG("Hit-cell distance in x-y is: " << proximity);

  // If the distance is positive then we are using the nearest cell rather than
  // are inside a cell If we are more than 2.25mm from the nearest cell we don't
  // create a hit to avoid the build-up of energy in edge cells For FCSV2
  // another hit can be created but with a cutoff to avoid looping, for
  // FastCaloGAN the rest of the hits in the layer will be scaled up by the
  // energy renormalization step.
  if (proximity < 2.25) {
    simulstate.deposit(cell.id(), hit.E());
  } else {
    hit.setXYZE(hit.x(), hit.y(), hit.z(), 0.0);
  }

  return FCSSuccess;
}

#pragma GCC diagnostic pop
