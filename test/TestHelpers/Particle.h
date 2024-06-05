#pragma once

#ifndef PARTICLE_H
#  define PARTICLE_H

#  include "CLHEP/Units/SystemOfUnits.h"
#  include "G4ThreeVector.hh"

namespace TestHelpers
{

struct Particle
{
  int pid;
  float ekin;
  G4ThreeVector pos;
  G4ThreeVector dir;

  friend void PrintTo(const Particle& particle, std::ostream* os)
  {
    // This method defines the name of the individual parametrized tests
    // testing::PrintToStringParamName() must be provided to
    // INSTANTIATE_TEST_SUITE_P as fourth argument

    // Eta as integer, e.g. 1.5 -> 150
    const int int_eta =
        static_cast<int>(std::round(particle.dir.eta() * 100.0));
    // Ekin as integer, e.g. 1.0 MeV -> 1
    const int int_ekin =
        static_cast<int>(std::round(particle.ekin * CLHEP::MeV));
    // Particle ID
    const int id = particle.pid;

    // Construct the string to be printed
    std::string print = "pid_" + std::to_string(id) + "_ekin_"
        + std::to_string(int_ekin) + "_eta_" + std::to_string(int_eta);

    // Replace '-' with 'm' in the string (e.g. -150 -> m150)
    std::replace(print.begin(), print.end(), '-', 'm');

    // Stream the string to the output stream
    *os << print;
  }
};

}  // namespace TestHelpers

#endif  // PARTICLE_H
