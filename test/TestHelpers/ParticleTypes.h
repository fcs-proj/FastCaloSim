#pragma once

#include <string_view>

namespace TestHelpers
{

struct ParticleType
{
  constexpr ParticleType(int pid, std::string_view label)
      : pid(pid)
      , label(label)
  {
  }

  const int pid;
  const std::string_view label;
};

constexpr auto make_particle(int pid, std::string_view label) -> ParticleType
{
  return ParticleType {pid, label};
}

namespace ParticleTypes
{
constexpr ParticleType Electron = make_particle(11, "e^-");
constexpr ParticleType Photon = make_particle(22, "\\gamma");
constexpr ParticleType Positron = make_particle(-11, "e^+");
constexpr ParticleType Proton = make_particle(2212, "p");
constexpr ParticleType Neutron = make_particle(2112, "n");
constexpr ParticleType PionPlus = make_particle(211, "\\pi^+");
constexpr ParticleType PionMinus = make_particle(-211, "\\pi^-");
}  // namespace ParticleTypes

}  // namespace TestHelpers
