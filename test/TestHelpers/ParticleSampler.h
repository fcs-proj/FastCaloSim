// Copyright (c) 2024 CERN for the benefit of the FastCaloSim project

#pragma once

#include <G4ThreeVector.hh>

#include "G4ThreeVector.hh"
#include "IDCaloBoundary.h"
#include "ParticleContainer.h"

namespace TestHelpers
{

struct IntersectionResult
{
  bool intersects;
  G4ThreeVector intersection;
};

class ParticleSampler
{
  // Pure vitural constructor
public:
  ParticleSampler() = default;
  virtual ~ParticleSampler() = default;

  virtual auto generate(int pid, float ekin, float eta) -> Particle = 0;
  virtual auto discreteEtaSample(int pid, float ekin, std::vector<float> etas)
      -> ParticleContainer = 0;
  virtual auto uniformEtaSample(int pid,
                                float ekin,
                                float min_eta,
                                float max_eta,
                                float step_size) -> ParticleContainer = 0;
};

class IDCaloBoundarySampler : public ParticleSampler
{
public:
  IDCaloBoundarySampler()
      : m_IDCaloBoundary()
  {
    // We set the thickness of the boundary to zero for the intersection
    // calculation
    m_IDCaloBoundary.setThickness(0);
  };

  inline auto generate(int pid, float ekin, float eta) -> Particle
  {
    // Detector assumed to be phi-symmetric, so pick random phi
    float phi = 0;
    float theta = 2 * std::atan(std::exp(-eta));
    G4ThreeVector direction;
    direction.setREtaPhi(1, eta, phi);
    direction = direction.unit();

    // Check for intersection of particle with ID-Calo boundary
    for (auto& cylinder : m_IDCaloBoundary.getCylinders()) {
      if (direction.eta() < 0) {
        // IDCaloBoundary cylinders are only defined for positive z
        // need to flip cylinders if eta < 0
        cylinder.second.zmin = -cylinder.second.zmax;
        cylinder.second.zmax = -cylinder.second.zmin;
      }

      auto result = intersects(cylinder.second, direction);
      if (result.intersects) {
        Particle particle =
            Particle {pid, ekin, result.intersection, direction};
        // Move particle slightly inside the boundary
        // This is currently a hack as particles do not seem
        // to be affected by magnetic field otherwise
        // Interestingly, magnetic field seems to get activated for eta<0 ?
        // TODO: discuss with Geant4 experts
        small_move_to_ID(particle, cylinder.second);
        return particle;
      }
    }

    // If no intersection was found, raise an exception
    throw std::runtime_error("No intersection found with ID-Calo boundary");
  };

  inline auto discreteEtaSample(int pid, float ekin, std::vector<float> etas)
      -> ParticleContainer
  {
    ParticleContainer particles;
    for (const auto& eta : etas) {
      Particle particle = generate(pid, ekin, eta);
      particles.add(particle);
    }
    return particles;
  };

  inline auto uniformEtaSample(int pid,
                               float ekin,
                               float min_eta,
                               float max_eta,
                               float step_size) -> ParticleContainer
  {
    if (min_eta > max_eta) {
      throw std::runtime_error("min_eta must be less than or equal to max_eta");
    }
    if (step_size < 0 || step_size > (max_eta - min_eta)) {
      throw std::runtime_error(
          "step_size must be zero or greater than zero and less than or equal "
          "to "
          "max_eta - min_eta");
    }

    // Calculate the number of steps
    // If min_eta == max_eta, we still want to sample at least once
    int n_steps = 1;
    if (!(min_eta == max_eta)) {
      n_steps =
          static_cast<int>(std::floor((max_eta - min_eta) / step_size)) + 1;
    }

    std::vector<float> etas;
    etas.reserve(n_steps);

    for (int i = 0; i < n_steps; ++i) {
      etas.emplace_back(min_eta + i * step_size);
    }
    return discreteEtaSample(pid, ekin, etas);
  };

private:
  IDCaloBoundary m_IDCaloBoundary;

  static inline auto intersects(Cylinder& cylinder, const G4ThreeVector& dir)
      -> IntersectionResult
  {
    IntersectionResult result;
    result.intersects = false;

    // Calculate intersection in z-direction for positive z
    double z = cylinder.zmax;
    double r = z * std::tan(dir.theta());

    if (r >= cylinder.rmin && r <= cylinder.rmax) {
      // Intersection within cylinder bounds
      double x = r * std::cos(dir.phi());
      double y = r * std::sin(dir.phi());

      result.intersects = true;
      result.intersection = G4ThreeVector(x, y, z);
      return result;
    }

    // Calculate intersection in z-direction for negative z
    z = cylinder.zmin;
    r = z * std::tan(dir.theta());

    if (r >= cylinder.rmin && r <= cylinder.rmax) {
      // Intersection within cylinder bounds
      double x = r * std::cos(dir.phi());
      double y = r * std::sin(dir.phi());

      result.intersects = true;
      result.intersection = G4ThreeVector(x, y, z);
      return result;
    }

    // Check intersection in r-direction (radial bounds)
    r = cylinder.rmax;
    z = r / std::tan(dir.theta());

    if (z >= cylinder.zmin && z <= cylinder.zmax) {
      double x = r * std::cos(dir.phi());
      double y = r * std::sin(dir.phi());

      result.intersects = true;
      result.intersection = G4ThreeVector(x, y, z);
      return result;
    }

    return result;
  }

  static inline auto small_move_to_ID(Particle& particle,
                                      Cylinder& cylinder,
                                      float epsilon = 0.01 * CLHEP::mm) -> void
  {
    if (cylinder.barrel) {
      float new_r = particle.pos.perp() - epsilon;
      particle.pos.setPerp(new_r);

    } else {
      if (particle.pos.eta() > 0) {
        float new_z = particle.pos.z() - epsilon;
        particle.pos.setZ(new_z);
      } else {
        float new_z = particle.pos.z() + epsilon;
        particle.pos.setZ(new_z);
      }
    }
  }
};

}  // namespace TestHelpers
