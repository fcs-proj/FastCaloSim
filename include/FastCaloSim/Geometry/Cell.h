#pragma once

#include <cassert>

#include <fmt/core.h>
#include <fmt/format.h>

#include "CLHEP/Vector/ThreeVector.h"

using Vector3D = CLHEP::Hep3Vector;

// Base class for cells
class Cell
{
private:
  /// @brief Unique identifier for the cell
  long long m_id;
  /// @brief Position of the cell in 3D space (x, y, z, eta, phi, r)
  Vector3D m_pos;
  /// @brief Layer of the cell
  long long m_layer;
  /// @brief Flag to indicate if the cell is in the barrel
  bool m_isBarrel;
  /// @brief Flags to indicate the coordinate system of the (cuboid) cell
  bool m_isXYZ, m_isEtaPhiR, m_isEtaPhiZ;
  /// @brief Cell sizes
  double m_dx, m_dy, m_dz, m_deta, m_dphi, m_dr;

public:
  Cell(long long id,
       Vector3D pos,
       long long layer,
       bool isBarrel,
       bool isXYZ,
       bool isEtaPhiR,
       bool isEtaPhiZ,
       double dx,
       double dy,
       double dz,
       double deta,
       double dphi,
       double dr)

      : m_id(id)
      , m_pos(std::move(pos))
      , m_layer(layer)
      , m_isBarrel(isBarrel)
      , m_isXYZ(isXYZ)
      , m_isEtaPhiR(isEtaPhiR)
      , m_isEtaPhiZ(isEtaPhiZ)
      , m_dx(dx)
      , m_dy(dy)
      , m_dz(dz)
      , m_deta(deta)
      , m_dphi(dphi)
      , m_dr(dr)
  {
  }

  virtual ~Cell() = default;

  // Direct accessors
  auto inline id() const -> long long { return m_id; }
  auto inline x() const -> double { return m_pos.x(); }
  auto inline y() const -> double { return m_pos.y(); }
  auto inline z() const -> double { return m_pos.z(); }
  auto inline eta() const -> double { return m_pos.eta(); }
  auto inline phi() const -> double { return m_pos.phi(); }
  auto inline r() const -> double { return m_pos.perp(); }
  auto inline layer() const -> long long { return m_layer; }
  auto inline isBarrel() const -> bool { return m_isBarrel; }

  auto inline isXYZ() const -> bool { return m_isXYZ; }
  auto inline isEtaPhiR() const -> bool { return m_isEtaPhiR; }
  auto inline isEtaPhiZ() const -> bool { return m_isEtaPhiZ; }

  auto inline dx() const -> double
  {
    assert(m_isXYZ && "Cell is not in XYZ coordinate system");
    return m_dx;
  }
  auto inline dy() const -> double
  {
    assert(m_isXYZ && "Cell is not in XYZ coordinate system");
    return m_dy;
  }

  auto inline dz() const -> double
  {
    assert(m_isXYZ
           || m_isEtaPhiZ && "Cell is not in XYZ or EtaPhiZ coordinate system");
    return m_dz;
  }

  auto inline deta() const -> double
  {
    assert(m_isEtaPhiR
           || m_isEtaPhiZ
               && "Cell is not in EtaPhiR or EtaPhiZ coordinate system");
    return m_deta;
  }

  auto inline dphi() const -> double
  {
    assert(m_isEtaPhiR
           || m_isEtaPhiZ
               && "Cell is not in EtaPhiR or EtaPhiZ coordinate system");
    return m_dphi;
  }

  auto inline dr() const -> double
  {
    assert(m_isEtaPhiR && "Cell is not in EtaPhiR coordinate system");
    return m_dr;
  }

  static auto norm_angle(double angle) -> double
  {
    angle = std::fmod(angle + M_PI, 2.0 * M_PI);
    if (angle < 0) {
      angle += 2.0 * M_PI;
    }
    return angle - M_PI;
  }

  /// @brief Calculates the 2D proximity between a hit and the boundary of a
  /// cell.
  ///
  /// This method returns the signed distance between the given hit and the cell
  /// boundary:
  /// - A negative value indicates that the hit is inside the cell, with the
  /// magnitude representing how deep within the cell the hit is.
  /// - A positive value indicates that the hit is outside the cell, with the
  /// magnitude representing the distance from the cell boundary.
  ///
  template<typename T>
  auto inline boundary_proximity(const T& hit) const -> double
  {
    if (m_isXYZ) {
      double delta_x = std::abs(hit.x() - m_pos.x());
      double delta_y = std::abs(hit.y() - m_pos.y());

      return std::max(delta_x - m_dx, delta_y - m_dy);
    }

    if (m_isEtaPhiR || m_isEtaPhiZ) {
      double delta_eta = std::abs(hit.eta() - m_pos.eta());
      double delta_phi = std::abs(norm_angle(hit.phi() - m_pos.phi()));

      return std::max(delta_eta - m_deta, delta_phi - m_dphi);
    }
  }

  template<typename T>
  auto inline is_inside(const T& hit) const -> bool
  {
    return boundary_proximity(hit) < 0;
  }

  // Overload the << operator to allow direct cell printout
  // with std::cout<<cell<<std::endl;
  friend auto operator<<(std::ostream& os, const Cell& cell) -> std::ostream&
  {
    os << fmt::format(
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
        cell.m_id,
        cell.m_layer,
        cell.m_pos.x(),
        cell.m_pos.y(),
        cell.m_pos.z(),
        cell.m_pos.phi(),
        cell.m_pos.eta(),
        cell.m_pos.perp());

    if (cell.m_isXYZ) {
      os << fmt::format(
          "Cell size:\n"
          "  dx: {:10.2f}\n"
          "  dy: {:10.2f}\n"
          "  dz: {:10.2f}\n"
          "---------------------------\n\n",
          cell.dx(),
          cell.dy(),
          cell.dz());
    }

    if (cell.m_isEtaPhiR) {
      os << fmt::format(
          "Cell size:\n"
          "  dphi: {:10.2f}\n"
          "  deta: {:10.2f}\n"
          "  dr: {:10.2f}\n"
          "---------------------------\n\n",
          cell.dphi(),
          cell.deta(),
          cell.dr());
    }

    if (cell.m_isEtaPhiZ) {
      os << fmt::format(
          "Cell size:\n"
          "  dphi: {:10.2f}\n"
          "  deta: {:10.2f}\n"
          "  dz: {:10.2f}\n"
          "---------------------------\n\n",
          cell.dphi(),
          cell.deta(),
          cell.dz());
    }
    return os;
  }
};
