// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#pragma once

#include <cstddef>

#include <G4ThreeVector.hh>

#include "FastCaloSim/Definitions/FloatCompare.h"
#include "IDCaloBoundary.h"
#include "ParticleContainer.h"

namespace TestHelpers
{

struct IntersectionResult
{
  bool intersects {false};
  G4ThreeVector intersection;
};

class ParticleSampler
{
  // Pure vitural constructor
public:
  ParticleSampler() = default;
  virtual ~ParticleSampler() = default;
  ParticleSampler(const ParticleSampler&) = default;
  ParticleSampler(ParticleSampler&&) = default;
  auto operator=(const ParticleSampler&) -> ParticleSampler& = default;
  auto operator=(ParticleSampler&&) -> ParticleSampler& = default;

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

  {
    // We set the thickness of the boundary to zero for the intersection
    // calculation
    m_IDCaloBoundary.setThickness(0);
  }

  auto generate(int pid, float ekin, float eta) -> Particle override
  {
    // Detector assumed to be phi-symmetric, so pick random phi
    float const phi = 0;
    G4ThreeVector direction;
    direction.setREtaPhi(1, static_cast<double>(eta), static_cast<double>(phi));
    direction = direction.unit();

    // Check for intersection of particle with ID-Calo boundary
    for (const auto& [name, boundary_cylinder] :
         m_IDCaloBoundary.getCylinders())
    {
      Cylinder cylinder = boundary_cylinder;
      if (direction.eta() < 0) {
        // IDCaloBoundary cylinders are only defined for positive z: mirror
        // this copy of the cylinder to negative z
        const double zmin = cylinder.zmin;
        cylinder.zmin = -cylinder.zmax;
        cylinder.zmax = -zmin;
      }

      auto result = intersects(cylinder, direction);
      if (result.intersects) {
        return Particle {.pid = pid,
                         .ekin = ekin,
                         .pos = result.intersection,
                         .dir = direction};
      }
    }

    // If no intersection was found, raise an exception
    throw std::runtime_error("No intersection found with ID-Calo boundary");
  }

  auto discreteEtaSample(int pid, float ekin, std::vector<float> etas)
      -> ParticleContainer override
  {
    ParticleContainer particles;
    for (const auto& eta : etas) {
      Particle const particle = generate(pid, ekin, eta);
      particles.add(particle);
    }
    return particles;
  }

  auto uniformEtaSample(int pid,
                        float ekin,
                        float min_eta,
                        float max_eta,
                        float step_size) -> ParticleContainer override
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
    if (!FastCaloSim::exactly_equal(min_eta, max_eta)) {
      n_steps =
          static_cast<int>(std::floor((max_eta - min_eta) / step_size)) + 1;
    }

    std::vector<float> etas;
    etas.reserve(static_cast<std::size_t>(n_steps));

    for (int i = 0; i < n_steps; ++i) {
      etas.emplace_back(min_eta + static_cast<float>(i) * step_size);
    }
    return discreteEtaSample(pid, ekin, etas);
  }

private:
  IDCaloBoundary m_IDCaloBoundary;

  static auto intersects(Cylinder& cylinder, const G4ThreeVector& dir)
      -> IntersectionResult
  {
    IntersectionResult result;
    result.intersects = false;

    // Calculate intersection in z-direction for positive z
    double z = cylinder.zmax;
    double r = z * std::tan(dir.theta());

    if (r >= cylinder.rmin && r <= cylinder.rmax) {
      // Intersection within cylinder bounds
      double const x = r * std::cos(dir.phi());
      double const y = r * std::sin(dir.phi());

      result.intersects = true;
      result.intersection = G4ThreeVector(x, y, z);
      return result;
    }

    // Calculate intersection in z-direction for negative z
    z = cylinder.zmin;
    r = z * std::tan(dir.theta());

    if (r >= cylinder.rmin && r <= cylinder.rmax) {
      // Intersection within cylinder bounds
      double const x = r * std::cos(dir.phi());
      double const y = r * std::sin(dir.phi());

      result.intersects = true;
      result.intersection = G4ThreeVector(x, y, z);
      return result;
    }

    // Check intersection in r-direction (radial bounds)
    r = cylinder.rmax;
    z = r / std::tan(dir.theta());

    if (z >= cylinder.zmin && z <= cylinder.zmax) {
      double const x = r * std::cos(dir.phi());
      double const y = r * std::sin(dir.phi());

      result.intersects = true;
      result.intersection = G4ThreeVector(x, y, z);
      return result;
    }

    return result;
  }
};

}  // namespace TestHelpers
