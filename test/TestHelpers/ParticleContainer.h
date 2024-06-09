#pragma once

#include <vector>

#include <G4ThreeVector.hh>

#include "Particle.h"

namespace TestHelpers
{

class ParticleContainer
{
public:
  // Default constructor
  ParticleContainer() = default;
  // Constructor with vector of particles
  explicit ParticleContainer(const std::vector<Particle>& particles)
      : m_particles(particles)
  {
  }

  // Add a particle to the container
  inline void addParticle(const Particle& particle)
  {
    m_particles.emplace_back(particle);
  }

  // Add a vector of particles to the container
  void addParticles(const std::vector<Particle>& particles)
  {
    m_particles.insert(m_particles.end(), particles.begin(), particles.end());
  }

  // Add a particle based on pid, eta, ekin, phi, pos
  inline void addParticle(int pid, float ekin, RPhiEta pos, RPhiEta dir)
  {
    Particle particle = {pid, ekin, pos.get(), dir.get().unit()};
    m_particles.emplace_back(particle);
  }

  inline void clear() { m_particles.clear(); }
  inline auto get() const -> const std::vector<Particle>&
  {
    return m_particles;
  }

private:
  std::vector<Particle> m_particles;
};

}  // namespace TestHelpers
