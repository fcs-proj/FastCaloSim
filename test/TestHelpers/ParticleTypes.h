// Copyright (c) 2024 CERN for the benefit of the FastCaloSim project

#pragma once

#include <string_view>

namespace TestHelpers
{

struct ParticleType
{
  constexpr ParticleType(int pid, std::string_view name, std::string_view label)
      : pid(pid)
      , name(name)
      , label(label)
  {
  }

  const int pid;
  const std::string_view name;
  const std::string_view label;
};

constexpr auto make_particle(int pid,
                             std::string_view name,
                             std::string_view label) -> ParticleType
{
  return ParticleType {pid, name, label};
}

namespace ParticleTypes
{
constexpr ParticleType Electron = make_particle(11, "electron", "e^-");
constexpr ParticleType Photon = make_particle(22, "photon", "\\gamma");
constexpr ParticleType Positron = make_particle(-11, "positron", "e^+");
constexpr ParticleType Proton = make_particle(2212, "proton", "p");
constexpr ParticleType Neutron = make_particle(2112, "neutron", "n");
constexpr ParticleType PionPlus = make_particle(211, "pi_plus", "\\pi^+");
constexpr ParticleType PionMinus = make_particle(-211, "pi_minus", "\\pi^-");
}  // namespace ParticleTypes

}  // namespace TestHelpers
