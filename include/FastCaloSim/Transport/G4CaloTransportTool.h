// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#pragma once

#include <memory>
#include <string>
#include <vector>

#include <FastCaloSim/FastCaloSim_export.h>

class G4Track;
class G4VPhysicalVolume;
class G4FieldTrack;
class G4Navigator;
class G4PropagatorInField;

/// @class G4CaloTransportTool
/// @brief A tool which transports particles through the Geant4 geometry.
///
/// The tool holds a Geant4 navigator and propagator, which keep the state of
/// the particle being transported. Like any Geant4 object with such state, it
/// belongs to one thread: use one instance per thread, constructed on the
/// thread that uses it, e.g. as a member of a per-thread fast simulation model.
/// The world volume it navigates is shared by all instances.
///
/// @author Joshua Falco Beirer <joshua.falco.beirer@cern.ch>
///
class FASTCALOSIM_EXPORT G4CaloTransportTool
{
public:
  /// @brief Transport particles through the given world volume, e.g. the one
  /// of createSimplifiedWorld() or the world of the full Geant4 geometry.
  explicit G4CaloTransportTool(G4VPhysicalVolume* world);
  ~G4CaloTransportTool();

  G4CaloTransportTool(const G4CaloTransportTool&) = delete;
  auto operator=(const G4CaloTransportTool&) -> G4CaloTransportTool& = delete;

  /// @brief Place the world volume of the simplified geometry, whose logical
  /// volumes must have been loaded (e.g. from GDML) into the Geant4 geometry
  /// stores. Call it once, on the master thread, and pass the result to the
  /// tools of all threads. The volume is owned by the Geant4 geometry store.
  static auto createSimplifiedWorld(const std::string& worldLogicalVolumeName)
      -> G4VPhysicalVolume*;

  /// @brief Transport the track through the geometry and return its steps
  auto transport(const G4Track& G4InputTrack) -> std::vector<G4FieldTrack>;

  /// @brief Set the name of the volume until which to transport the particle
  void setTransportLimitVolume(const std::string& transportLimitVolume)
  {
    m_transportLimitVolume = transportLimitVolume;
  }

  /// @brief Set the maximum number of steps allowed in particle transport
  void setMaxSteps(unsigned int maxSteps) { m_maxSteps = maxSteps; }

private:
  /// @brief Advance the track by one Geant4 step, starting in the volume the
  /// navigator located the track in last.
  void doStep(G4FieldTrack& fieldTrack, G4VPhysicalVolume* currentPhysVol);

  /// Navigator of this tool in the shared world volume
  std::unique_ptr<G4Navigator> m_navigator;

  /// Propagator of charged particles in the magnetic field, using m_navigator
  std::unique_ptr<G4PropagatorInField> m_propagator;

  /// Name of the volume in which to stop tracking the particle
  std::string m_transportLimitVolume = "Envelope";

  /// Maximum number of steps allowed in particle transport
  unsigned int m_maxSteps = 100;
};
