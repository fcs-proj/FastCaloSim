// Copyright (c) 2024 CERN for the benefit of the FastCaloSim project

#ifndef G4ATLASTOOLS_G4CALOTRANSPORTTOOL_H
#define G4ATLASTOOLS_G4CALOTRANSPORTTOOL_H

#include <FastCaloSim/FastCaloSim_export.h>

#include "FastCaloSim/Transport/ThreadLocalHolder.h"
#include "G4PropagatorInField.hh"

class G4Track;
class G4VPhysicalVolume;
class G4FieldTrack;

/// @class G4CaloTransportTool
/// @brief A tool which transports particles through the Geant4 geometry.
///
/// @author Joshua Falco Beirer <joshua.falco.beirer@cern.ch>
///
class FASTCALOSIM_EXPORT G4CaloTransportTool
{
public:
  G4CaloTransportTool();
  ~G4CaloTransportTool();

  // Initialize propagator for the current thread
  void initializePropagator();

  // Transport input track through the geometry
  auto transport(const G4Track& G4InputTrack) -> std::vector<G4FieldTrack>;

  //--------------------------------------------------------------------------
  //  Setter methods for configuration
  //--------------------------------------------------------------------------
  /// @brief Turn on/off the use of simplified geometry.
  void setUseSimplifiedGeo(bool useSimplifiedGeo)
  {
    m_useSimplifiedGeo = useSimplifiedGeo;
  }

  /// @brief Set the name of the logical volume for the simplified world
  void setSimplifiedWorldLogName(const std::string& simplifiedWorldLogName)
  {
    m_simplifiedWorldLogName = simplifiedWorldLogName;
  }

  /// @brief Set the name of the volume until which to transport the particle
  void setTransportLimitVolume(const std::string& transportLimitVolume)
  {
    m_transportLimitVolume = transportLimitVolume;
  }

  /// @brief Set the maximum number of steps allowed in particle transport
  void setMaxSteps(unsigned int maxSteps) { m_maxSteps = maxSteps; }

private:
  /// @brief Obtain the world volume in which particle transport is performed
  auto getWorldVolume() -> G4VPhysicalVolume*;

  /// @brief Create and return a new G4PropagatorInField.
  auto makePropagator() -> G4PropagatorInField*;

  /// @brief Advance the track by one Geant4 step in the geometry
  void doStep(G4FieldTrack& fieldTrack);

  /// Pointer to the physical volume of the world (either simplified or full
  /// geometry)
  G4VPhysicalVolume* m_worldVolume {};

  /// Whether to use simplified geometry for particle transport
  bool m_useSimplifiedGeo = true;

  /// Name of the logical volume of the simplified world (as defined in GDML)
  std::string m_simplifiedWorldLogName = "WorldLog";

  /// Name of the volume in which to stop tracking the particle
  std::string m_transportLimitVolume = "Envelope";

  /// Maximum number of steps allowed in particle transport
  unsigned int m_maxSteps = 100;

  /// Thread-local holder for G4PropagatorInField instances
  thread_utils::ThreadLocalHolder<G4PropagatorInField> m_propagatorHolder;
};

#endif  // G4ATLASTOOLS_G4CALOTRANSPORTTOOL_H
