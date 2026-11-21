/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

//
// BFieldH8Grid.cxx
//
// Masahiro Morii, Harvard University
//
#include <array>
#include <cstddef>
#include <string>

#include "BFieldH8Grid.h"

BFieldH8Grid::BFieldH8Grid() = default;

void BFieldH8Grid::readMap(std::istream& input)
{
  // read the magnet header line
  std::string name;
  double cor {};
  m_n[0] = 0;
  input >> name;
  for (int& i : m_n) {
    input >> i;
  }
  for (std::size_t i = 0; i < 3; i++) {
    input >> m_min[i] >> m_max[i];
  }
  for (double& i : m_d) {
    input >> i;
  }
  input >> cor;
  // anything read?
  if (m_n[0] == 0) {
    return;
  }
  // convert unit to mm
  const double meter(1000.);  // m in mm
  for (std::size_t i = 0; i < 3; i++) {
    m_min[i] *= meter;
    m_max[i] *= meter;
    m_d[i] *= meter;
  }
  // prepare space for the field data
  int const nxyz = m_n[0] * m_n[1] * m_n[2];
  for (auto& i : m_B) {
    i.resize(static_cast<std::size_t>(nxyz));
  }
  // read the field data
  std::array<double, 3> xyz {};
  std::array<double, 3> B {};
  for (int k = 0; k < nxyz; k++) {
    input >> xyz[0] >> xyz[1] >> xyz[2] >> B[0] >> B[1] >> B[2];
    // adjust unit and scale
    const double cm(10.);  // cm in mm
    const double tesla(0.001);  // T in kT
    for (std::size_t i = 0; i < 3; i++) {
      xyz[i] *= cm;
      B[i] *= -cor * tesla;
    }
    // find the grid number
    std::array<int, 3> j {};
    for (std::size_t i = 0; i < 3; i++) {
      const double grid_pos =
          (xyz[i] - m_min[i]) / (m_max[i] - m_min[i]) * (m_n[i] - 1);
      // NOLINTNEXTLINE(bugprone-incorrect-roundings): as in Athena
      j[i] = static_cast<int>(grid_pos + 0.5);
    }
    // store
    const int index = j[0] + m_n[0] * (j[1] + m_n[1] * j[2]);
    const auto ixyz = static_cast<std::size_t>(index);
    for (std::size_t i = 0; i < 3; i++) {
      m_B[i][ixyz] = B[i];
    }
  }
  // apply offset to the range
  for (std::size_t i = 0; i < 3; i++) {
    m_min[i] += m_d[i];
    m_max[i] += m_d[i];
  }
}

void BFieldH8Grid::getB(const double* xyz, double* B, double* deriv) const
{
  if (!defined() || !inside(xyz)) {
    B[0] = B[1] = B[2] = 0.0;
    if (deriv != nullptr) {
      for (int j = 0; j < 9; j++) {
        deriv[j] = 0.0;
      }
    }
    return;
  }
  // find the grid index
  std::array<double, 3> invunit {};
  std::array<int, 3> j {};
  std::array<double, 3> f {};
  std::array<double, 3> g {};
  for (std::size_t i = 0; i < 3; i++) {
    invunit[i] = (m_n[i] - 1) / (m_max[i] - m_min[i]);
    double const a = (xyz[i] - m_min[i]) * invunit[i];
    j[i] = static_cast<int>(a);
    f[i] = a - j[i];
    g[i] = 1.0 - f[i];
  }
  // indices of the 8 corners
  const auto nx = static_cast<std::size_t>(m_n[0]);
  const int n_xy = m_n[0] * m_n[1];
  const auto nxy = static_cast<std::size_t>(n_xy);
  const int index = j[0] + m_n[0] * (j[1] + m_n[1] * j[2]);
  const auto ixyz0 = static_cast<std::size_t>(index);
  const std::array<std::size_t, 8> ixyz {ixyz0,
                                         ixyz0 + 1,
                                         ixyz0 + nx,
                                         ixyz0 + nx + 1,
                                         ixyz0 + nxy,
                                         ixyz0 + nxy + 1,
                                         ixyz0 + nxy + nx,
                                         ixyz0 + nxy + nx + 1};
  // interpolate field values
  for (std::size_t i = 0; i < 3; i++) {  // Bx, By, Bz
    B[i] = g[2]
            * (g[1] * (g[0] * m_B[i][ixyz[0]] + f[0] * m_B[i][ixyz[1]])
               + f[1] * (g[0] * m_B[i][ixyz[2]] + f[0] * m_B[i][ixyz[3]]))
        + f[2]
            * (g[1] * (g[0] * m_B[i][ixyz[4]] + f[0] * m_B[i][ixyz[5]])
               + f[1] * (g[0] * m_B[i][ixyz[6]] + f[0] * m_B[i][ixyz[7]]));
  }
  // derivatives
  if (deriv != nullptr) {
    for (std::size_t i = 0; i < 3; i++) {  // Bx, By, Bz
      deriv[i * 3] = (g[2]
                          * (g[1] * (m_B[i][ixyz[1]] - m_B[i][ixyz[0]])
                             + f[1] * (m_B[i][ixyz[3]] - m_B[i][ixyz[2]]))
                      + f[2]
                          * (g[1] * (m_B[i][ixyz[5]] - m_B[i][ixyz[4]])
                             + f[1] * (m_B[i][ixyz[7]] - m_B[i][ixyz[6]])))
          * invunit[0];
      deriv[i * 3 + 1] = (g[2]
                              * (g[0] * (m_B[i][ixyz[2]] - m_B[i][ixyz[0]])
                                 + f[0] * (m_B[i][ixyz[3]] - m_B[i][ixyz[1]]))
                          + f[2]
                              * (g[0] * (m_B[i][ixyz[6]] - m_B[i][ixyz[4]])
                                 + f[0] * (m_B[i][ixyz[7]] - m_B[i][ixyz[5]])))
          * invunit[1];
      deriv[i * 3 + 2] = (g[1]
                              * (g[0] * (m_B[i][ixyz[4]] - m_B[i][ixyz[0]])
                                 + f[0] * (m_B[i][ixyz[5]] - m_B[i][ixyz[1]]))
                          + f[1]
                              * (g[0] * (m_B[i][ixyz[6]] - m_B[i][ixyz[2]])
                                 + f[0] * (m_B[i][ixyz[7]] - m_B[i][ixyz[3]])))
          * invunit[2];
    }
  }
}

void BFieldH8Grid::setOffset(const double* dxyz)
{
  // update the range by the change in offset
  for (std::size_t i = 0; i < 3; i++) {
    m_min[i] += dxyz[i] - m_d[i];
    m_max[i] += dxyz[i] - m_d[i];
    m_d[i] = dxyz[i];
  }
}
