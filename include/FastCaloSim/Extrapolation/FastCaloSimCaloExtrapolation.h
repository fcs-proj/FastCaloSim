// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#pragma once

#include <vector>

#include <FastCaloSim/FastCaloSim_export.h>

#include "FastCaloSim/Core/MLogging.h"

class CaloGeo;
class G4FieldTrack;
class TFCSExtrapolationState;
class TFCSTruthState;

/// @class FastCaloSimCaloExtrapolation
/// @brief Extrapolates a particle, given by the steps of its transport, to the
/// ID-calo boundary and to the entrance, middle and exit of each calorimeter
/// layer.
///
/// Each of these surfaces is approximated by a cylinder around the beam axis.
/// The particle is extrapolated to it along the straight lines between its
/// transport steps.
class FASTCALOSIM_EXPORT FastCaloSimCaloExtrapolation
    : public FastCaloSim::MLogging
{
public:
  /// @brief Set the geometry of the calorimeter layers
  void set_geometry(CaloGeo* geo) { m_geo = geo; }

  /// @brief Fill result with the positions of the particle on the ID-calo
  /// boundary and in the calorimeter layers. The layers are left empty if the
  /// particle does not reach the calorimeter.
  void extrapolate(TFCSExtrapolationState& result,
                   const TFCSTruthState* truth,
                   const std::vector<G4FieldTrack>& caloSteps) const;

private:
  CaloGeo* m_geo = nullptr;
};
