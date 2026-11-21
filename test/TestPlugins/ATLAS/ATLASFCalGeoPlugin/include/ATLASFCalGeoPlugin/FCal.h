// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#pragma once

#include <algorithm>
#include <array>
#include <cfloat>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>

#include <RtypesCore.h>
#include <TString.h>

#include "FCAL_ChannelMap.h"
#include "FastCaloSim/Definitions/FloatCompare.h"
#include "FastCaloSim/Geometry/CaloGeo.h"

class FCal : public CaloGeo
{
public:
  /// @brief Default constructor
  FCal() = default;

  // Virtual destructor to ensure proper cleanup
  ~FCal() override = default;

  /// @brief Set the geo pointer
  auto set_geo(CaloGeo* geo) { m_geo = geo; }

  /// @brief Get the cell at a specific position
  auto get_cell(std::size_t layer, const Position& pos) const
      -> const Cell& override
  {
    unsigned long long const cell_id = this->get_fcal_cell_id(
        static_cast<int>(layer), pos.x(), pos.y(), pos.z());

    if (cell_id == static_cast<unsigned long long>(-1)) {
      static Cell const invalid_cell;
      return invalid_cell;
    }
    // This will get the cell directly from the cell store
    return m_geo->get_cell(cell_id);
  }

  /// @brief Load the real FCal geometry from a file
  auto load(const std::array<std::string, 3>& fileNames) -> void
  {
    std::vector<std::unique_ptr<std::istream>> electrodes;
    std::vector<std::unique_ptr<std::streambuf>> electrodesBuf;
    electrodes.reserve(3);
    electrodesBuf.reserve(3);

    for (uint16_t i = 0; i < 3; i++) {
      const std::string& file = fileNames.at(i);

      std::unique_ptr<std::ifstream> directStream =
          std::make_unique<std::ifstream>(file);
      if (!directStream->is_open()) {
        std::cout << "Error: Could not open " << file << '\n';
        throw std::runtime_error("Could not open file");
      }
      electrodes.push_back(std::move(directStream));
    }

    int thisTubeId;
    int thisTubeI;
    int thisTubeJ;
    // int    thisTubeID;
    // int    thisTubeMod;
    double thisTubeX;
    double thisTubeY;
    TString tubeName;

    // int second_column;
    std::string seventh_column;
    std::string eight_column;
    int ninth_column;

    for (int imodule = 1; imodule <= 3; imodule++) {
      // std::cout << "Loading FCal electrode #" << imodule << std::endl;

      std::istream& electrode =
          *electrodes[static_cast<std::size_t>(imodule - 1)];
      while (1) {
        electrode >> tubeName;
        if (electrode.eof()) {
          break;
        }
        electrode >> thisTubeId;  // ?????
        electrode >> thisTubeI;
        electrode >> thisTubeJ;
        electrode >> thisTubeX;
        electrode >> thisTubeY;
        electrode >> seventh_column;
        electrode >> eight_column;
        electrode >> ninth_column;

        tubeName.ReplaceAll("'", "");
        std::string const tubeNamestring = tubeName.Data();

        std::istringstream tileStream1(std::string(tubeNamestring, 1, 1));
        std::istringstream tileStream2(std::string(tubeNamestring, 3, 2));
        std::istringstream tileStream3(std::string(tubeNamestring, 6, 3));
        int a1 = 0;
        int a2 = 0;
        int a3 = 0;
        if (tileStream1) {
          tileStream1 >> a1;
        }
        if (tileStream2) {
          tileStream2 >> a2;
        }
        if (tileStream3) {
          tileStream3 >> a3;
        }

        std::stringstream const s;

        m_FCal_ChannelMap.add_tube(tubeNamestring,
                                   imodule,
                                   thisTubeId,
                                   thisTubeI,
                                   thisTubeJ,
                                   thisTubeX,
                                   thisTubeY,
                                   seventh_column);
      }
    }

    m_FCal_ChannelMap.finish();  // Creates maps

    electrodes.clear();
    electrodesBuf.clear();

    this->calculateFCalRminRmax();
  }

  /// @brief Calculates the minimum and maximum radii of the FCal cells
  void calculateFCalRminRmax()
  {
    m_FCal_rmin.resize(3, static_cast<double>(FLT_MAX));
    m_FCal_rmax.resize(3, 0.);

    double x(0.);
    double y(0.);
    double r(0.);
    for (int imap = 1; imap <= 3; imap++) {
      for (auto it = m_FCal_ChannelMap.begin(imap);
           it != m_FCal_ChannelMap.end(imap);
           it++)
      {
        const auto idx = static_cast<std::size_t>(imap - 1);
        x = static_cast<double>(it->second.x());
        y = static_cast<double>(it->second.y());
        r = std::sqrt(x * x + y * y);
        m_FCal_rmin[idx] = std::min(r, m_FCal_rmin[idx]);
        m_FCal_rmax[idx] = std::max(r, m_FCal_rmax[idx]);
      }
    }
  }

  auto getClosestFCalCellIndex(
      int layer, float x, float y, int& ieta, int& iphi, int* steps) const
      -> bool
  {
    double const rmin = m_FCal_rmin[static_cast<std::size_t>(layer - 21)];
    double const rmax = m_FCal_rmax[static_cast<std::size_t>(layer - 21)];
    int const isam = layer - 20;
    double a = 1.;
    const double b = 0.01;
    const int nmax = 100;
    int i = 0;

    // x * x + y * y is evaluated in float, the square root in double
    const double r = sqrt(static_cast<double>(x * x + y * y));
    if (FastCaloSim::exactly_equal(r, 0.)) {
      return false;
    }
    const double r_inverse = 1. / r;

    // Scale a float coordinate in double precision, round back to float
    auto scaled = [](double s, float v) -> float
    { return static_cast<float>(s * static_cast<double>(v)); };

    if ((r / rmax) > (rmin * r_inverse)) {
      x = static_cast<float>(static_cast<double>(x) * rmax * r_inverse);
      y = static_cast<float>(static_cast<double>(y) * rmax * r_inverse);
      while ((!m_FCal_ChannelMap.getTileID(
                 isam, scaled(a, x), scaled(a, y), ieta, iphi))
             && i < nmax)
      {
        a -= b;
        i++;
      }
    } else {
      x = static_cast<float>(static_cast<double>(x) * rmin * r_inverse);
      y = static_cast<float>(static_cast<double>(y) * rmin * r_inverse);
      while ((!m_FCal_ChannelMap.getTileID(
                 isam, scaled(a, x), scaled(a, y), ieta, iphi))
             && i < nmax)
      {
        a += b;
        i++;
      }
    }
    if (steps) {
      *steps = i + 1;
    }
    return i < nmax;
  }

  /// @brief Returns the best matching FCAL cell for a hit in a specific layer
  auto get_fcal_cell_id(int layer,
                        double x,
                        double y,
                        double z,
                        int* steps = nullptr) const -> unsigned long long
  {
    // The FCal modules are the layers 21, 22 and 23
    int const isam = layer - 20;
    if (isam < 1 || isam > 3) {
      std::cout << "Error: Layer " << layer << " is not an FCal layer!" << '\n';
      return static_cast<unsigned long long>(-1);
    }
    int iphi(-100000);
    int ieta(-100000);
    Long64_t const mask1[] {0x34, 0x34, 0x35};
    Long64_t const mask2[] {0x36, 0x36, 0x37};
    // The channel map works in single precision
    const auto xf = static_cast<float>(x);
    const auto yf = static_cast<float>(y);
    bool found = m_FCal_ChannelMap.getTileID(isam, xf, yf, ieta, iphi);
    if (steps && found) {
      *steps = 0;
    }
    if (!found) {
      // Hit is not matched with any FCal cell! Looking for
      // the closest cell
      found = getClosestFCalCellIndex(layer, xf, yf, ieta, iphi, steps);
    }
    if (!found) {
      std::cout << "Error: Unable to find the closest FCal cell!" << '\n';
      return static_cast<unsigned long long>(-1);
    }
    Long64_t id = (ieta << 5) + 2 * iphi;
    if (isam == 2) {
      id += (8 << 8);
    }

    if (z > 0) {
      id += (mask2[isam - 1] << 12);
    } else {
      id += (mask1[isam - 1] << 12);
    }

    id = id << 44;

    return static_cast<unsigned long long>(id);
  }

protected:
  // The constructor argument is set according to
  // https://atlas-sw-doxygen.web.cern.ch/atlas-sw-doxygen/atlas_main--Doxygen/docs/html/d6/d72/FCALChannelMapBuilder_8cxx_source.html
  FCAL_ChannelMap m_FCal_ChannelMap {0};
  std::vector<double> m_FCal_rmin, m_FCal_rmax;

  // Pointer to geo as we use the FastCaloSim cell store
  CaloGeo* m_geo;
};
