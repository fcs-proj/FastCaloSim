// Copyright (c) 2024 CERN for the benefit of the FastCaloSim project

/*
  Access point for particle data (currently only electric charge) needed during
  simulation.

  FastCaloSim ships a small built-in table that covers the particles produced by
  the standard parametrisations. Experiments that need fuller particle coverage
  can provide their own implementation by deriving from ParticleData::Provider
  and installing it once at start-up via ParticleData::setProvider, e.g.
  wrapping a generic source such as the AtlasPID / TruthUtils helpers:
  https://gitlab.cern.ch/atlas/athena/-/blob/main/Generators/TruthUtils/TruthUtils/AtlasPID.h?ref_type=heads

  All call sites use the free ParticleData::charge() function, which dispatches
  to the installed provider (or the built-in default if none was installed).
*/

#pragma once

#include <memory>

#include <FastCaloSim/FastCaloSim_export.h>

namespace ParticleData
{

/// Interface for supplying particle data to FastCaloSim.
///
/// Experiments implement this to plug in their own particle database. The
/// interface is intentionally minimal; further per-particle quantities (mass,
/// name, stability, ...) can be added here as the need arises, and every
/// implementation gets them through the single injection point.
class FASTCALOSIM_EXPORT Provider
{
public:
  virtual ~Provider() = default;

  /// Electric charge of the particle with the given PDG id, in units of e.
  /// Implementations should throw if the pdgID is not supported.
  virtual auto charge(int pdgID) const -> double = 0;
};

/// Install the provider used by all subsequent particle data queries.
///
/// Intended to be called once during experiment initialization, before any
/// simulation runs. Passing nullptr restores the built-in default provider.
/// Not safe to call concurrently with simulation.
FASTCALOSIM_EXPORT auto setProvider(std::shared_ptr<const Provider> provider)
    -> void;

/// The currently installed provider, or the built-in default if none was set.
/// Returned by shared_ptr so the provider stays alive for the caller even if
/// setProvider replaces it concurrently.
FASTCALOSIM_EXPORT auto getProvider() -> std::shared_ptr<const Provider>;

/// Convenience accessor: electric charge of the given PDG id, in units of e.
inline auto charge(int pdgID) -> double
{
  return getProvider()->charge(pdgID);
}

}  // namespace ParticleData
