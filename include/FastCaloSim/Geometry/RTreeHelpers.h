// Copyright (c) 2025 CERN for the benefit of the FastCaloSim project
#pragma once

#include <array>
#include <cmath>
#include <vector>

#include "FastCaloSim/Geometry/Cell.h"

namespace RTreeHelpers
{
enum class CoordinateSystem
{
  XYZ,
  EtaPhiR,
  EtaPhiZ,
  RPhiZ,
  Undefined
};

// Normalize phi angle to [-π, π)
inline double norm_phi(double angle)
{
  angle = std::fmod(angle + M_PI, 2.0 * M_PI);
  if (angle < 0)
    angle += 2.0 * M_PI;
  return angle - M_PI;
}

// Create bounding boxes for eta-phi coordinates
inline std::vector<std::array<double, 4>> build_eta_phi_boxes(double eta,
                                                              double phi,
                                                              double deta,
                                                              double dphi)
{
  double half_deta = deta / 2.0;
  double half_dphi = dphi / 2.0;
  double etamin = eta - half_deta;
  double etamax = eta + half_deta;
  double phimin = norm_phi(phi - half_dphi);
  double phimax = norm_phi(phi + half_dphi);

  std::vector<std::array<double, 4>> boxes;
  if (phimin > phimax) {
    // Handle boxes that cross the -π/π boundary
    boxes.push_back(std::array<double, 4> {etamin, -M_PI, etamax, phimax});
    boxes.push_back(std::array<double, 4> {etamin, phimin, etamax, M_PI});
  } else {
    boxes.push_back(std::array<double, 4> {etamin, phimin, etamax, phimax});
  }
  return boxes;
}

// Create bounding box for x-y coordinates
inline std::array<double, 4> build_xy_box(double x,
                                          double y,
                                          double dx,
                                          double dy)
{
  double half_dx = dx / 2.0;
  double half_dy = dy / 2.0;
  return std::array<double, 4> {
      x - half_dx, y - half_dy, x + half_dx, y + half_dy};
}
}  // namespace RTreeHelpers
