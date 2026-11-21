// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#pragma once

#include <map>
#include <string>

#include <CLHEP/Units/SystemOfUnits.h>

namespace TestHelpers
{

struct Cylinder
{
  double rmin;
  double rmax;
  double zmin;
  double zmax;
  bool barrel;
};

class IDCaloBoundary
{
public:
  IDCaloBoundary() { recalculate(); }

  void setThickness(float thickness)
  {
    cylinderThickness = thickness;
    recalculate();
  }

  auto getThickness() const -> float { return cylinderThickness; }

  auto getCylinders() -> std::map<std::string, Cylinder>&
  {
    return m_cylinders;
  }

private:
  float cylinderThickness {10 * CLHEP::mm};

  std::map<std::string, Cylinder> m_cylinders;

  void recalculate()
  {
    // The cylinder dimensions are stored in double precision
    const auto thickness = static_cast<double>(cylinderThickness);
    m_cylinders["barrel"] = {.rmin = 1148 * CLHEP::mm,
                             .rmax = 1148 * CLHEP::mm + thickness,
                             .zmin = 0.0 * CLHEP::mm,
                             .zmax = 3550.0 * CLHEP::mm,
                             .barrel = true};

    m_cylinders["barrelEndcap"] = {
        .rmin = 120 * CLHEP::mm,
        .rmax = m_cylinders["barrel"].rmax,
        .zmin = m_cylinders["barrel"].zmax,
        .zmax = m_cylinders["barrel"].zmax + thickness,
        .barrel = false};

    m_cylinders["innerBeamPipe"] = {
        .rmin = m_cylinders["barrelEndcap"].rmin,
        .rmax = m_cylinders["barrelEndcap"].rmin + thickness,
        .zmin = m_cylinders["barrelEndcap"].zmax,
        .zmax = 4587.0 * CLHEP::mm,
        .barrel = true};

    m_cylinders["innerBeamPipeEndcap"] = {
        .rmin = 41 * CLHEP::mm,
        .rmax = m_cylinders["innerBeamPipe"].rmin,
        .zmin = m_cylinders["innerBeamPipe"].zmax,
        .zmax = m_cylinders["innerBeamPipe"].zmax + thickness,
        .barrel = false};

    m_cylinders["outerBeamPipe"] = {
        .rmin = m_cylinders["innerBeamPipeEndcap"].rmin,
        .rmax = m_cylinders["innerBeamPipeEndcap"].rmin + thickness,
        .zmin = m_cylinders["innerBeamPipeEndcap"].zmax,
        .zmax = 6783 * CLHEP::mm,
        .barrel = true};
  }
};

}  // namespace TestHelpers
