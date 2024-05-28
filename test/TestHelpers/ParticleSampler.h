#ifndef PARTICLE_SAMPLER_H
#define PARTICLE_SAMPLER_H

#include <vector>

#include <G4ThreeVector.hh>

#include "Particle.h"

namespace TestHelpers::ParticleSampler
{

/// @brief Generate a list of particles from lists of pids and etas.
/// @param pid_list A list of particle IDs you want to generate.
/// @param eta_list A list of eta values you want to generate.
/// @param ekin The kinetic energy of the particles.
/// @param phi The phi value of the particles.
/// @param pos The primary vertex position of the particles.
/// @return A vector of generated particles.
inline auto generate_particles(const std::vector<int>& pid_list,
                               const std::vector<float>& eta_list,
                               float ekin,
                               float phi,
                               const G4ThreeVector& pos)
    -> std::vector<TestHelpers::Particle>
{
  std::vector<TestHelpers::Particle> particles;
  for (const int pid : pid_list) {
    for (const float eta : eta_list) {
      G4ThreeVector dir;
      dir.setRhoPhiEta(1.0, phi, eta);
      TestHelpers::Particle particle = {pid, ekin, pos, dir};
      particles.emplace_back(particle);
    }
  }
  return particles;
}

}  // namespace TestHelpers::ParticleSampler

#endif  // PARTICLE_SAMPLER_H
