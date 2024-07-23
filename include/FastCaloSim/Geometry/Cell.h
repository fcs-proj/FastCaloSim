#pragma once

#include <string>

#include <fmt/core.h>
#include <fmt/format.h>
#include <pv/polymorphic_variant.hpp>

// Base class for cells
class CellBase
{
private:
  long long m_id;
  double m_x, m_y, m_z, m_phi, m_eta, m_r;
  long long m_layer;
  bool m_isXYZ, m_isEtaPhiR, m_isEtaPhiZ;

public:
  CellBase(long long id,
           double x,
           double y,
           double z,
           double phi,
           double eta,
           double r,
           int layer,
           bool isXYZ,
           bool isEtaPhiR,
           bool isEtaPhiZ)
      : m_id(id)
      , m_x(x)
      , m_y(y)
      , m_z(z)
      , m_phi(phi)
      , m_eta(eta)
      , m_r(r)
      , m_layer(layer)
      , m_isXYZ(isXYZ)
      , m_isEtaPhiR(isEtaPhiR)
      , m_isEtaPhiZ(isEtaPhiZ)
  {
  }

  // Direct accessors
  auto inline id() const -> long long { return m_id; }
  auto inline x() const -> double { return m_x; }
  auto inline y() const -> double { return m_y; }
  auto inline z() const -> double { return m_z; }
  auto inline phi() const -> double { return m_phi; }
  auto inline eta() const -> double { return m_eta; }
  auto inline r() const -> double { return m_r; }
  auto inline layer() const -> long long { return m_layer; }
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

  virtual void print() const
  {
    std::string output = fmt::format(
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
        m_id,
        m_layer,
        m_x,
        m_y,
        m_z,
        m_phi,
        m_eta,
        m_r);
    fmt::print("{}", output);
  }
};

// Derived class for Cartesian cells (dx, dy, dz)
class XYZCell : public CellBase
{
private:
  double m_dx, m_dy, m_dz;

public:
  XYZCell(long long id,
          double x,
          double y,
          double z,
          double phi,
          double eta,
          double r,
          long long layer,
          double dx,
          double dy,
          double dz)
      : CellBase(id, x, y, z, phi, eta, r, layer, true, false, false)
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

  void print() const override
  {
    CellBase::print();
    std::string output = fmt::format(
        "Cell size:\n"
        "  dx: {:10.2f}\n"
        "  dy: {:10.2f}\n"
        "  dz: {:10.2f}\n"
        "---------------------------\n\n",
        m_dx,
        m_dy,
        m_dz);
    fmt::print("{}", output);
  }
};

// Derived class for Cylindrical cells (deta, dphi, dr)
class EtaPhiRCell : public CellBase
{
private:
  double m_deta, m_dphi, m_dr;

public:
  EtaPhiRCell(long long id,
              double x,
              double y,
              double z,
              double phi,
              double eta,
              double r,
              long long layer,
              double deta,
              double dphi,
              double dr)
      : CellBase(id, x, y, z, phi, eta, r, layer, false, true, false)
      , m_deta(deta)
      , m_dphi(dphi)
      , m_dr(dr)
  {
  }

  auto inline deta() const -> double override { return m_deta; }
  auto inline dphi() const -> double { return m_dphi; }
  auto inline dr() const -> double { return m_dr; }

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

  void print() const override
  {
    CellBase::print();
    std::string output = fmt::format(
        "Cell size:\n"
        "  dphi: {:10.2f}\n"
        "  deta: {:10.2f}\n"
        "  dr: {:10.2f}\n"
        "---------------------------\n\n",
        m_dphi,
        m_deta,
        m_dr);
    fmt::print("{}", output);
  }
};

// Derived class for Cylindrical cells with z dimension (deta, dphi, dz)
class EtaPhiZCell : public CellBase
{
private:
  double m_deta, m_dphi, m_dz;

public:
  EtaPhiZCell(long long id,
              double x,
              double y,
              double z,
              double phi,
              double eta,
              double r,
              long long layer,
              double deta,
              double dphi,
              double dz)
      : CellBase(id, x, y, z, phi, eta, r, layer, false, false, true)
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

  void print() const override
  {
    CellBase::print();
    std::string output = fmt::format(
        "Cell size:\n"
        "  dphi: {:10.2f}\n"
        "  deta: {:10.2f}\n"
        "  dz: {:10.2f}\n"
        "---------------------------\n\n",
        m_dphi,
        m_deta,
        m_dz);
    fmt::print("{}", output);
  }
};

// A cell is any variant of XYZCell, EtaPhiRCell, EtaPhiZCell
// Instead of std::variant, we use a polymorphic variant to get easy access of
// base class methods for variants See
// https://github.com/Krzmbrzl/polymorphic_variant for more information
using Cell =
    pv::polymorphic_variant<CellBase, XYZCell, EtaPhiRCell, EtaPhiZCell>;
