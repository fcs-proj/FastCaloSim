// Copyright (c) 2024 CERN for the benefit of the FastCaloSim project

#include <cstdlib>
#include <memory>
#include <stdexcept>
#include <string>

#include "FastCaloSim/Definitions/ParticleData.h"

namespace ParticleData
{

namespace
{

/// Built-in particle data used when no experiment provider is installed.
/// Covers only the particles produced by the standard parametrisations.
class DefaultProvider final : public Provider
{
public:
  auto charge(int pdgID) const -> double override
  {
    if (pdgID == 11 || pdgID == 211 || pdgID == 2212)
      return 1.;
    if (pdgID == -11 || pdgID == -211 || pdgID == -2212)
      return -1.;
    if (pdgID == 22 || std::abs(pdgID) == 2112 || pdgID == 111)
      return 0.;

    throw std::runtime_error("ParticleData: unsupported pdgID: "
                             + std::to_string(pdgID));
  }
};

/// The installed provider, lazily defaulting to DefaultProvider. Defined in a
/// single translation unit and reached through the exported accessors below so
/// that there is exactly one instance across the shared-library boundary.
auto activeProvider() -> std::shared_ptr<const Provider>&
{
  static std::shared_ptr<const Provider> provider =
      std::make_shared<const DefaultProvider>();
  return provider;
}

}  // namespace

auto setProvider(std::shared_ptr<const Provider> provider) -> void
{
  activeProvider() = provider ? std::move(provider)
                              : std::make_shared<const DefaultProvider>();
}

auto getProvider() -> std::shared_ptr<const Provider>
{
  return activeProvider();
}

}  // namespace ParticleData
