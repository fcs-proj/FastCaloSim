#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-equal"

/*
  Copyright (C) 2002-2018 CERN for the benefit of the ATLAS collaboration
*/

#include "FastCaloSim/Core/TFCSHitCellMappingFCal.h"

#include "FastCaloSim/Core/TFCSSimulationState.h"
#include "FastCaloSim/Geometry/ICaloGeometry.h"

//=============================================
//======= TFCSHitCellMappingFCal =========
//=============================================

FCSReturnCode TFCSHitCellMappingFCal::simulate_hit(
    Hit& hit,
    TFCSSimulationState& simulstate,
    const TFCSTruthState* /*truth*/,
    const TFCSExtrapolationState* /*extrapol*/)
{
  int cs = calosample();
  float distance;
  const CaloDetDescrElement* cellele =
      m_geo->getFCalDDE(cs, hit.x(), hit.y(), hit.z(), &distance);
  ATH_MSG_DEBUG("HIT: cellele=" << cellele << " E=" << hit.E() << " cs=" << cs
                                << " x=" << hit.x() << " y=" << hit.y()
                                << " z=" << hit.z());

  /// protection against cases where hits cannot be matched to a FCal cell
  if ((hit.x() == 0 && hit.y() == 0) || cellele == nullptr) {
    ATH_MSG_WARNING(
        "TFCSLateralShapeParametrizationHitCellMapping::simulate_hit: cellele="
        << cellele << " E=" << hit.E() << " cs=" << cs << " eta=" << hit.eta()
        << " phi=" << hit.phi());
    return (FCSReturnCode)(FCSRetry + 5);  // retry simulation up to 5 times
  }

  // If the distance is positive then we are using the nearest cell rather than
  // are inside a cell If we are more than 2.25mm from the nearest cell we don't
  // create a hit to avoid the build-up of energy in edge cells For FCSV2
  // another hit can be created but with a cutoff to avoid looping, for
  // FastCaloGAN the rest of the hits in the layer will be scaled up by the
  // energy renormalization step.
  if (distance < 2.25) {
    simulstate.deposit(cellele, hit.E());
  } else {
    hit.setXYZE(hit.x(), hit.y(), hit.z(), 0.0);
  }

  return FCSSuccess;
}

#pragma GCC diagnostic pop
