#pragma once

#include <fmt/core.h>
#include <fmt/format.h>
#include <pv/polymorphic_variant.hpp>

#include "CLHEP/Vector/ThreeVector.h"

using Vector3D = CLHEP::Hep3Vector;

// Base class for cells
class CellBase
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

public:
  CellBase(long long id,
           Vector3D pos,
           int layer,
           bool isBarrel,
           bool isXYZ,
           bool isEtaPhiR,
           bool isEtaPhiZ)
      : m_id(id)
      , m_pos(std::move(pos))
      , m_layer(layer)
      , m_isBarrel(isBarrel)
      , m_isXYZ(isXYZ)
      , m_isEtaPhiR(isEtaPhiR)
      , m_isEtaPhiZ(isEtaPhiZ)
  {
  }

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

  // Virtual accessors to be implemented by the daughter classes
  virtual auto dx() const -> double = 0;
  virtual auto dy() const -> double = 0;
  virtual auto dz() const -> double = 0;
  virtual auto deta() const -> double = 0;
  virtual auto dphi() const -> double = 0;
  virtual auto dr() const -> double = 0;

  virtual ~CellBase() = default;

  // Overload the << operator to allow direct cell printout
  // with std::cout<<cell<<std::endl;
  friend auto operator<<(std::ostream& os, const CellBase& cell)
      -> std::ostream&
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
      return os;
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
      return os;
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
      return os;
    }
  }
};

// Derived class for Cartesian cells (dx, dy, dz)
class XYZCell : public CellBase
{
private:
  double m_dx, m_dy, m_dz;

public:
  XYZCell(long long id,
          Vector3D pos,
          long long layer,
          bool isBarrel,
          double dx,
          double dy,
          double dz)
      : CellBase(id, pos, layer, isBarrel, true, false, false)
      , m_dx(dx)
      , m_dy(dy)
      , m_dz(dz)
  {
  }

  auto inline dx() const -> double override { return m_dx; }
  auto inline dy() const -> double override { return m_dy; }
  auto inline dz() const -> double override { return m_dz; }

  auto inline deta() const -> double override
  {
    throw std::runtime_error("Cannot call deta on XYZ cell");
  }
  auto inline dphi() const -> double override
  {
    throw std::runtime_error("Cannot call dphi on XYZ cell");
  }
  auto inline dr() const -> double override
  {
    throw std::runtime_error("Cannot call dr on XYZ cell");
  }
};

// Derived class for Cylindrical cells (deta, dphi, dr)
class EtaPhiRCell : public CellBase
{
private:
  double m_deta, m_dphi, m_dr;

public:
  EtaPhiRCell(long long id,
              Vector3D pos,
              long long layer,
              bool isBarrel,
              double deta,
              double dphi,
              double dr)
      : CellBase(id, pos, layer, isBarrel, false, true, false)
      , m_deta(deta)
      , m_dphi(dphi)
      , m_dr(dr)
  {
  }

  auto inline deta() const -> double override { return m_deta; }
  auto inline dphi() const -> double override { return m_dphi; }
  auto inline dr() const -> double override { return m_dr; }

  auto inline dx() const -> double override
  {
    throw std::runtime_error("Cannot call dx on EtaPhiR cell");
  }
  auto inline dy() const -> double override
  {
    throw std::runtime_error("Cannot call dy on EtaPhiR cell");
  }
  auto inline dz() const -> double override
  {
    throw std::runtime_error("Cannot call dz on EtaPhiR cell");
  }
};

// Derived class for Cylindrical cells with z dimension (deta, dphi, dz)
class EtaPhiZCell : public CellBase
{
private:
  double m_deta, m_dphi, m_dz;

public:
  EtaPhiZCell(long long id,
              Vector3D pos,
              long long layer,
              bool isBarrel,
              double deta,
              double dphi,
              double dz)
      : CellBase(id, pos, layer, isBarrel, false, false, true)
      , m_deta(deta)
      , m_dphi(dphi)
      , m_dz(dz)
  {
  }

  auto inline deta() const -> double override { return m_deta; }
  auto inline dphi() const -> double override { return m_dphi; }
  auto inline dz() const -> double override { return m_dz; }

  auto inline dx() const -> double override
  {
    throw std::runtime_error("Cannot call dx on EtaPhiZ cell");
  }
  auto inline dy() const -> double override
  {
    throw std::runtime_error("Cannot call dy on EtaPhiZ cell");
  }
  auto inline dr() const -> double override
  {
    throw std::runtime_error("Cannot call dr on EtaPhiZ cell");
  }
};

// A cell is any variant of XYZCell, EtaPhiRCell, EtaPhiZCell
// Instead of std::variant, we use a polymorphic variant to get easy access of
// base class methods for variants See
// https://github.com/Krzmbrzl/polymorphic_variant for more information
using Cell =
    pv::polymorphic_variant<CellBase, XYZCell, EtaPhiRCell, EtaPhiZCell>;
