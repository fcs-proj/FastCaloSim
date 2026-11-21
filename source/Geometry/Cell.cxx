// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#include <format>
#include <ostream>

#include "FastCaloSim/Geometry/Cell.h"

auto operator<<(std::ostream& os, const Cell& cell) -> std::ostream&
{
  os << std::format(
      "Cell details:\n"
      "---------------------------\n"
      "ID: {}\n"
      "Layer: {}\n"
      "Coordinates:\n"
      "  x: {:10.2f}\n"
      "  y: {:10.2f}\n"
      "  z: {:10.2f}\n"
      "  phi: {:10.2f}\n"
      "  eta: {:10.2f}\n"
      "  r: {:10.2f}\n",
      cell.id(),
      cell.layer(),
      cell.x(),
      cell.y(),
      cell.z(),
      cell.phi(),
      cell.eta(),
      cell.r());

  if (cell.isXYZ()) {
    os << std::format(
        "Cell size:\n"
        "  dx: {:10.2f}\n"
        "  dy: {:10.2f}\n"
        "  dz: {:10.2f}\n"
        "---------------------------\n\n",
        cell.dx(),
        cell.dy(),
        cell.dz());
  }
  if (cell.isEtaPhiR()) {
    os << std::format(
        "Cell size:\n"
        "  dphi: {:10.2f}\n"
        "  deta: {:10.2f}\n"
        "  dr: {:10.2f}\n"
        "---------------------------\n\n",
        cell.dphi(),
        cell.deta(),
        cell.dr());
  }
  if (cell.isEtaPhiZ()) {
    os << std::format(
        "Cell size:\n"
        "  dphi: {:10.2f}\n"
        "  deta: {:10.2f}\n"
        "  dz: {:10.2f}\n"
        "---------------------------\n\n",
        cell.dphi(),
        cell.deta(),
        cell.dz());
  }
  if (cell.isRPhiZ()) {
    os << std::format(
        "Cell size:\n"
        "  dphi: {:10.2f}\n"
        "  dz: {:10.2f}\n"
        "  dr: {:10.2f}\n"
        "---------------------------\n\n",
        cell.dphi(),
        cell.dz(),
        cell.dr());
  }
  return os;
}
