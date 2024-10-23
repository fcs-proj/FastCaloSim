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
  std::vector<G4FieldTrack> transport(const G4Track& G4InputTrack);

private:
  // Get the world volume in which the particle transport is performed
  G4VPhysicalVolume* getWorldVolume();
  // Create and return a new propagator
  G4PropagatorInField* makePropagator();
  // Advance track by single Geant4 step in geometry
  void doStep(G4FieldTrack& fieldTrack);
  // Pointer to the physical volume of the world (either simplified or full
  // geometry)
  G4VPhysicalVolume* m_worldVolume {};

  // Whether to use simplified geometry for particle transport
  bool m_useSimplifiedGeo = true;
  // Name of the logical volume of the simplified world as defined in the loaded
  // GDML file
  std::string m_simplifiedWorldLogName = "WorldLog";
  // Name of volume until which the particle is tracked in transport
  std::string m_transportLimitVolume = "Envelope";
  // Maximum number of steps in particle transport
  unsigned int m_maxSteps = 100;
  // Thread local holder for propagators
  thread_utils::ThreadLocalHolder<G4PropagatorInField> m_propagatorHolder;

};  // class G4CaloTransportTool

#endif  // G4ATLASTOOLS_G4CALOTRANSPORTTOOL_H
