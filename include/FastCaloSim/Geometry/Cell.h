#pragma once

#include <cassert>

#include <fmt/core.h>
#include <fmt/format.h>

struct Position
{
  double m_x, m_y, m_z;
  double m_eta, m_phi, m_r;
  auto x() const -> double { return m_x; }
  auto y() const -> double { return m_y; }
  auto z() const -> double { return m_z; }
  auto eta() const -> double { return m_eta; }
  auto phi() const -> double { return m_phi; }
  auto r() const -> double { return m_r; }
};

// Base class for cells
class Cell
{
private:
  /// @brief Unique identifier for the cell
  long long m_id;
  /// @brief Position of the cell in 3D space (x, y, z, eta, phi, r)
  Position m_pos;
  /// @brief Layer of the cell
  long long m_layer;
  /// @brief Flag to indicate if the cell is in the barrel
  bool m_isBarrel;
  /// @brief Flags to indicate the coordinate system of the (cuboid) cell
  bool m_isXYZ, m_isEtaPhiR, m_isEtaPhiZ;
  /// @brief Cell sizes
  double m_dx, m_dy, m_dz, m_deta, m_dphi, m_dr;
  /// @brief Flag to invalidate cell
  bool m_is_valid = true;

public:
  enum SubPos
  {
    MID = 0,  // middle
    ENT = 1,  // entrance
    EXT = 2  // exit
  };

  Cell(long long id,
       Position pos,
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

  // Allow easy creation of an invalid cell
  Cell()
      : m_id(-1)
      , m_pos(Position {0, 0, 0, 0, 0, 0})
      , m_layer(-1)
      , m_isBarrel(false)
      , m_isXYZ(false)
      , m_isEtaPhiR(false)
      , m_isEtaPhiZ(false)
      , m_dx(0)
      , m_dy(0)
      , m_dz(0)
      , m_deta(0)
      , m_dphi(0)
      , m_dr(0)
      , m_is_valid(false)
  {
  }

  virtual ~Cell() = default;

  // Cell validity
  auto inline is_valid() const -> bool { return m_is_valid; }
  auto inline invalidate() -> void { m_is_valid = false; }

  // Direct accessors
  auto inline id() const -> long long { return m_id; }
  auto inline x() const -> double { return m_pos.x(); }
  auto inline y() const -> double { return m_pos.y(); }
  auto inline z() const -> double { return m_pos.z(); }
  auto inline eta() const -> double { return m_pos.eta(); }
  auto inline phi() const -> double { return m_pos.phi(); }
  auto inline r() const -> double { return m_pos.r(); }
  auto inline layer() const -> long long { return m_layer; }
  auto inline isBarrel() const -> bool { return m_isBarrel; }

  auto inline isXYZ() const -> bool { return m_isXYZ; }
  auto inline isEtaPhiR() const -> bool { return m_isEtaPhiR; }
  auto inline isEtaPhiZ() const -> bool { return m_isEtaPhiZ; }

  /// @brief Cell size accessors
  /// Note: For the hit->cell matching, only the cell sizes
  /// for the coordinate system of the cell are used
  /// (e.g. for x,y,z cells, only dx, dy, dz are needed)
  /// However: for the extrapolation we currently also take into account
  /// the half width of the cell in the other coordinate systems
  /// e.g. for TileBar0 (eta,phi,r) we also use dz which is stored in the ROOT
  /// file In essence: if you call e.g. zent() on any cell it will add whatever
  /// dz is stored We might want to revisit the logic in the future
  /// to prevent misuse for now, we assert that stored cell widths are non-zero
  /// in the case e.g. zent() or zext() is called
  auto inline dx() const -> double { return m_dx; }
  auto inline dy() const -> double { return m_dy; }

  auto inline dz() const -> double { return m_dz; }

  auto inline deta() const -> double { return m_deta; }

  auto inline dphi() const -> double { return m_dphi; }

  auto inline dr() const -> double { return m_dr; }

  auto inline rent() const -> double
  {
    assert(m_dr > 0 && "rent() called on cell with dr <= 0. The half-width of the cell seems undefined.");
    return r() - m_dr * 0.5;
  }

  auto inline rext() const -> double
  {
    assert(m_dr > 0 && "rext() called on cell with dr <= 0.  The half-width of the cell seems undefined.");
    return r() + m_dr * 0.5;
  }

  auto inline zent() const -> double
  {
    assert(m_dz > 0 && "zent() called on cell with dz <= 0. The half-width of the cell seems undefined.");
    return z() < 0 ? z() + m_dz * 0.5 : z() - m_dz * 0.5;
  }

  auto inline zext() const -> double
  {
    assert(m_dz > 0 && "zext() called on cell with dz <= 0.  The half-width of the cell seems undefined.");
    return z() < 0 ? z() - m_dz * 0.5 : z() + m_dz * 0.5;
  }

  // only makes ense for barrel
  auto inline r(SubPos subpos) const -> double
  {
    switch (subpos) {
      case SubPos::ENT:
        return rent();
      case SubPos::MID:
        return r();
      case SubPos::EXT:
        return rext();
    }
  }

  auto inline z(SubPos subpos) const -> double
  {
    switch (subpos) {
      case SubPos::ENT:
        return zent();
      case SubPos::MID:
        return z();
      case SubPos::EXT:
        return zext();
    }
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
  /// This method returns the signed distance between the given hit and the
  /// cell boundary:
  /// - A negative value indicates that the hit is inside the cell, with the
  /// magnitude representing how deep within the cell the hit is.
  /// - A positive value indicates that the hit is outside the cell, with the
  /// magnitude representing the distance from the cell boundary.
  ///
  auto inline boundary_proximity(const Position& pos) const -> double
  {
    if (m_isXYZ) {
      double delta_x = std::abs(pos.x() - m_pos.x());
      double delta_y = std::abs(pos.y() - m_pos.y());

      return std::max(delta_x - m_dx, delta_y - m_dy);
    }

    if (m_isEtaPhiR || m_isEtaPhiZ) {
      double delta_eta = std::abs(pos.eta() - m_pos.eta());
      double delta_phi = std::abs(norm_angle(pos.phi() - m_pos.phi()));

      return std::max(delta_eta - m_deta, delta_phi - m_dphi);
    }

    std::runtime_error("Cell is not in a valid coordinate system");
    return 0;
  }

  auto inline is_inside(const Position& pos) const -> bool
  {
    return boundary_proximity(pos) < 0;
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
        cell.m_pos.r());

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
