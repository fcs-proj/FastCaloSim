#pragma once

#include <fmt/core.h>
#include <fmt/format.h>
#include <pv/polymorphic_variant.hpp>

// General struct to hold a position in the calorimeter
struct CaloPos
{
  double m_x, m_y, m_z, m_eta, m_phi, m_r;

  auto x() const -> double { return m_x; }
  auto y() const -> double { return m_y; }
  auto z() const -> double { return m_z; }
  auto eta() const -> double { return m_eta; }
  auto phi() const -> double { return m_phi; }
  auto r() const -> double { return m_r; }

  auto set_x(double x) -> void { m_x = x; }
  auto set_y(double y) -> void { m_y = y; }
  auto set_z(double z) -> void { m_z = z; }
  auto set_eta(double eta) -> void { m_eta = eta; }
  auto set_phi(double phi) -> void { m_phi = phi; }
  auto set_r(double r) -> void { m_r = r; }
};

// Base class for cells
class CellBase
{
private:
  long long m_id;
  CaloPos m_pos;
  long long m_layer;
  bool m_isBarrel;
  bool m_isXYZ, m_isEtaPhiR, m_isEtaPhiZ;

public:
  CellBase(long long id,
           CaloPos pos,
           int layer,
           bool isBarrel,
           bool isXYZ,
           bool isEtaPhiR,
           bool isEtaPhiZ)
      : m_id(id)
      , m_pos(pos)
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
  auto inline r() const -> double { return m_pos.r(); }
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
        cell.m_pos.r());
    return os;
  }
};

// Derived class for Cartesian cells (dx, dy, dz)
class XYZCell : public CellBase
{
private:
  double m_dx, m_dy, m_dz;

public:
  XYZCell(long long id,
          CaloPos pos,
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

  friend auto operator<<(std::ostream& os, const XYZCell& cell) -> std::ostream&
  {
    os << static_cast<const CellBase&>(cell);
    os << fmt::format(
        "Cell size:\n"
        "  dx: {:10.2f}\n"
        "  dy: {:10.2f}\n"
        "  dz: {:10.2f}\n"
        "---------------------------\n\n",
        cell.m_dx,
        cell.m_dy,
        cell.m_dz);
    return os;
  }
};

// Derived class for Cylindrical cells (deta, dphi, dr)
class EtaPhiRCell : public CellBase
{
private:
  double m_deta, m_dphi, m_dr;

public:
  EtaPhiRCell(long long id,
              CaloPos pos,
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

  friend auto operator<<(std::ostream& os, const EtaPhiRCell& cell)
      -> std::ostream&
  {
    os << static_cast<const CellBase&>(cell);
    os << fmt::format(
        "Cell size:\n"
        "  dphi: {:10.2f}\n"
        "  deta: {:10.2f}\n"
        "  dr: {:10.2f}\n"
        "---------------------------\n\n",
        cell.m_dphi,
        cell.m_deta,
        cell.m_dr);
    return os;
  }
};

// Derived class for Cylindrical cells with z dimension (deta, dphi, dz)
class EtaPhiZCell : public CellBase
{
private:
  double m_deta, m_dphi, m_dz;

public:
  EtaPhiZCell(long long id,
              CaloPos pos,
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

  friend auto operator<<(std::ostream& os, const EtaPhiZCell& cell)
      -> std::ostream&
  {
    os << static_cast<const CellBase&>(cell);
    os << fmt::format(
        "Cell size:\n"
        "  dphi: {:10.2f}\n"
        "  deta: {:10.2f}\n"
        "  dr: {:10.2f}\n"
        "---------------------------\n\n",
        cell.m_dphi,
        cell.m_deta,
        cell.m_dz);
    return os;
  }
};

// A cell is any variant of XYZCell, EtaPhiRCell, EtaPhiZCell
// Instead of std::variant, we use a polymorphic variant to get easy access of
// base class methods for variants See
// https://github.com/Krzmbrzl/polymorphic_variant for more information
using Cell =
    pv::polymorphic_variant<CellBase, XYZCell, EtaPhiRCell, EtaPhiZCell>;
