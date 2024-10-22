// Copyright (c) 2024 CERN for the benefit of the FastCaloSim project

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
  IDCaloBoundary()
      : cylinderThickness(10 * CLHEP::mm)
  {
    recalculate();
  }

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
  float cylinderThickness;

  std::map<std::string, Cylinder> m_cylinders;

  void recalculate()
  {
    m_cylinders["barrel"] = {1148 * CLHEP::mm,
                             1148 * CLHEP::mm + cylinderThickness,
                             0.0 * CLHEP::mm,
                             3550.0 * CLHEP::mm,
                             true};

    m_cylinders["barrelEndcap"] = {
        120 * CLHEP::mm,
        m_cylinders["barrel"].rmax,
        m_cylinders["barrel"].zmax,
        m_cylinders["barrel"].zmax + cylinderThickness,
        false};

    m_cylinders["innerBeamPipe"] = {
        m_cylinders["barrelEndcap"].rmin,
        m_cylinders["barrelEndcap"].rmin + cylinderThickness,
        m_cylinders["barrelEndcap"].zmax,
        4587.0 * CLHEP::mm,
        true};

    m_cylinders["innerBeamPipeEndcap"] = {
        41 * CLHEP::mm,
        m_cylinders["innerBeamPipe"].rmin,
        m_cylinders["innerBeamPipe"].zmax,
        m_cylinders["innerBeamPipe"].zmax + cylinderThickness,
        false};

    m_cylinders["outerBeamPipe"] = {
        m_cylinders["innerBeamPipeEndcap"].rmin,
        m_cylinders["innerBeamPipeEndcap"].rmin + cylinderThickness,
        m_cylinders["innerBeamPipeEndcap"].zmax,
        6783 * CLHEP::mm,
        true};
  }
};

}  // namespace TestHelpers
