// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#pragma once

#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <ostream>
#include <utility>

struct Position
{
  double m_x, m_y, m_z;
  double m_eta, m_phi, m_r;
  constexpr auto x() const -> double { return m_x; }
  constexpr auto y() const -> double { return m_y; }
  constexpr auto z() const -> double { return m_z; }
  constexpr auto eta() const -> double { return m_eta; }
  constexpr auto phi() const -> double { return m_phi; }
  constexpr auto r() const -> double { return m_r; }
};

struct CellData
{
  /// @brief Unique identifier for the cell
  uint64_t m_id;
  /// @brief Position of the cell in 3D space (x, y, z, eta, phi, r)
  Position m_pos;
  /// @brief Layer of the cell
  uint32_t m_layer;
  /// @brief Flag to indicate if the cell is in the barrel
  bool m_isBarrel;
  /// @brief Flags to indicate the coordinate system of the (cuboid) cell
  bool m_isXYZ, m_isEtaPhiR, m_isEtaPhiZ, m_isRPhiZ;
  /// @brief Cell sizes
  double m_dx, m_dy, m_dz, m_deta, m_dphi, m_dr;
  /// @brief Flag to invalidate cell
  bool m_is_valid;
};

static_assert(std::is_trivially_copyable<CellData>::value,
              "CellData must be trivially copiable");

class Cell
{
public:
  enum SubPos
  {
    MID = 0,  // middle
    ENT = 1,  // entrance
    EXT = 2  // exit
  };

  Cell()
  {
    m_data.m_id = static_cast<uint64_t>(-1);
    m_data.m_pos = {0, 0, 0, 0, 0, 0};
    m_data.m_layer = static_cast<uint32_t>(-1);
    m_data.m_isBarrel = false;
    m_data.m_isXYZ = m_data.m_isEtaPhiR = m_data.m_isEtaPhiZ =
        m_data.m_isRPhiZ = false;
    m_data.m_dx = m_data.m_dy = m_data.m_dz = 0;
    m_data.m_deta = m_data.m_dphi = m_data.m_dr = 0;
    m_data.m_is_valid = false;
  }

  Cell(uint64_t id,
       Position pos,
       uint32_t layer,
       bool isBarrel,
       bool isXYZ,
       bool isEtaPhiR,
       bool isEtaPhiZ,
       bool isRPhiZ,
       double dx,
       double dy,
       double dz,
       double deta,
       double dphi,
       double dr)
  {
    m_data.m_id = id;
    m_data.m_pos = pos;
    m_data.m_layer = layer;
    m_data.m_isBarrel = isBarrel;
    m_data.m_isXYZ = isXYZ;
    m_data.m_isEtaPhiR = isEtaPhiR;
    m_data.m_isEtaPhiZ = isEtaPhiZ;
    m_data.m_isRPhiZ = isRPhiZ;
    m_data.m_dx = dx;
    m_data.m_dy = dy;
    m_data.m_dz = dz;
    m_data.m_deta = deta;
    m_data.m_dphi = dphi;
    m_data.m_dr = dr;
    m_data.m_is_valid = true;
  }

  auto id() const -> uint64_t { return m_data.m_id; }
  auto x() const -> double { return m_data.m_pos.m_x; }
  auto y() const -> double { return m_data.m_pos.m_y; }
  auto z() const -> double { return m_data.m_pos.m_z; }
  auto eta() const -> double { return m_data.m_pos.m_eta; }
  auto phi() const -> double { return m_data.m_pos.m_phi; }
  auto r() const -> double { return m_data.m_pos.m_r; }
  auto layer() const -> uint32_t { return m_data.m_layer; }
  auto isBarrel() const -> bool { return m_data.m_isBarrel; }
  auto isXYZ() const -> bool { return m_data.m_isXYZ; }
  auto isEtaPhiR() const -> bool { return m_data.m_isEtaPhiR; }
  auto isEtaPhiZ() const -> bool { return m_data.m_isEtaPhiZ; }
  auto isRPhiZ() const -> bool { return m_data.m_isRPhiZ; }
  auto dx() const -> double { return m_data.m_dx; }
  auto dy() const -> double { return m_data.m_dy; }
  auto dz() const -> double { return m_data.m_dz; }
  auto deta() const -> double { return m_data.m_deta; }
  auto dphi() const -> double { return m_data.m_dphi; }
  auto dr() const -> double { return m_data.m_dr; }
  auto is_valid() const -> bool { return m_data.m_is_valid; }
  void invalidate() { m_data.m_is_valid = false; }
  // Raw access to the cell data
  auto raw() const -> const CellData& { return m_data; }
  auto raw() -> CellData& { return m_data; }

  auto rent() const -> double
  {
    assert(m_data.m_dr > 0 && "rent() called on cell with dr <= 0. The half-width of the cell seems undefined.");
    return r() - m_data.m_dr * 0.5;
  }
  auto rext() const -> double
  {
    assert(m_data.m_dr > 0 && "rext() called on cell with dr <= 0.  The half-width of the cell seems undefined.");
    return r() + m_data.m_dr * 0.5;
  }
  auto zent() const -> double
  {
    assert(m_data.m_dz > 0 && "zent() called on cell with dz <= 0. The half-width of the cell seems undefined.");
    return z() < 0 ? z() + m_data.m_dz * 0.5 : z() - m_data.m_dz * 0.5;
  }
  auto zext() const -> double
  {
    assert(m_data.m_dz > 0 && "zext() called on cell with dz <= 0.  The half-width of the cell seems undefined.");
    return z() < 0 ? z() - m_data.m_dz * 0.5 : z() + m_data.m_dz * 0.5;
  }

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
    std::unreachable();
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
    std::unreachable();
  }

  static auto norm_angle(double angle) -> double
  {
    constexpr double two_pi = 2.0 * M_PI;
    angle += M_PI;
    // Equivalent to angle = std::fmod(angle, two_pi), which is slow. For the
    // differences of two azimuths it removes at most one period, and that
    // subtraction is exact (Sterbenz lemma), so the result is identical.
    if (angle >= two_pi && angle < 2 * two_pi) {
      angle -= two_pi;
    } else if (!(angle > -two_pi && angle < two_pi)) {
      angle = std::fmod(angle, two_pi);
    }
    if (angle < 0) {
      angle += two_pi;
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
  auto boundary_proximity(const Position& pos) const -> double
  {
    if (m_data.m_isXYZ) {
      double dx = std::abs(pos.m_x - x());
      double dy = std::abs(pos.m_y - y());
      return std::max(dx - m_data.m_dx / 2, dy - m_data.m_dy / 2);
    }
    if (m_data.m_isEtaPhiR || m_data.m_isEtaPhiZ || m_data.m_isRPhiZ) {
      double d_eta = std::abs(pos.m_eta - eta());
      double d_phi = std::abs(norm_angle(pos.m_phi - phi()));
      return std::max(d_eta / m_data.m_deta, d_phi / m_data.m_dphi) - 0.5;
    }
    throw std::runtime_error("Invalid coordinate system");
  }

  auto is_inside(const Position& pos) const -> bool
  {
    return boundary_proximity(pos) < 0;
  }

  // Overload the << operator to allow direct cell printout
  // with std::cout<<cell<<std::endl;
  // Defined out of line so that the formatting code is not parsed by rootcling.
  friend auto operator<<(std::ostream& os, const Cell& cell) -> std::ostream&;

private:
  CellData m_data;
};
