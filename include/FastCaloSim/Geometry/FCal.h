#pragma once

#include <array>
#include <cfloat>
#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>

#include <RtypesCore.h>
#include <TString.h>

#include "FastCaloSim/Geometry/CaloGeo.h"
#include "FastCaloSim/Geometry/FCAL_ChannelMap.h"

class FCAL : public CaloGeo
{
public:
  /// @brief Default constructor
  FCAL() = default;

  // Virtual destructor to ensure proper cleanup
  virtual ~FCAL() = default;

  /// @brief Get the cell at a specific position
  auto get_cell_id(int layer, const Position& pos) const -> long long override
  {
    return this->get_fcal_cell_id(layer, pos.x(), pos.y(), pos.z());
  }

  /// @brief Load the real FCal geometry from a file
  auto load(const std::array<std::string, 3>& fileNames) -> bool
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
        std::cout << "Error: Could not open " << file << std::endl;
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

    int i;
    for (int imodule = 1; imodule <= 3; imodule++) {
      // std::cout << "Loading FCal electrode #" << imodule << std::endl;

      i = 0;
      while (1) {
        (*electrodes[imodule - 1]) >> tubeName;
        if (electrodes[imodule - 1]->eof())
          break;
        (*electrodes[imodule - 1]) >> thisTubeId;  // ?????
        (*electrodes[imodule - 1]) >> thisTubeI;
        (*electrodes[imodule - 1]) >> thisTubeJ;
        (*electrodes[imodule - 1]) >> thisTubeX;
        (*electrodes[imodule - 1]) >> thisTubeY;
        (*electrodes[imodule - 1]) >> seventh_column;
        (*electrodes[imodule - 1]) >> eight_column;
        (*electrodes[imodule - 1]) >> ninth_column;

        tubeName.ReplaceAll("'", "");
        std::string tubeNamestring = tubeName.Data();

        std::istringstream tileStream1(std::string(tubeNamestring, 1, 1));
        std::istringstream tileStream2(std::string(tubeNamestring, 3, 2));
        std::istringstream tileStream3(std::string(tubeNamestring, 6, 3));
        int a1 = 0, a2 = 0, a3 = 0;
        if (tileStream1)
          tileStream1 >> a1;
        if (tileStream2)
          tileStream2 >> a2;
        if (tileStream3)
          tileStream3 >> a3;

        std::stringstream s;

        m_FCal_ChannelMap.add_tube(tubeNamestring,
                                   imodule,
                                   thisTubeId,
                                   thisTubeI,
                                   thisTubeJ,
                                   thisTubeX,
                                   thisTubeY,
                                   seventh_column);

        i++;
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
    m_FCal_rmin.resize(3, FLT_MAX);
    m_FCal_rmax.resize(3, 0.);

    double x(0.), y(0.), r(0.);
    for (int imap = 1; imap <= 3; imap++)
      for (auto it = m_FCal_ChannelMap.begin(imap);
           it != m_FCal_ChannelMap.end(imap);
           it++)
      {
        x = it->second.x();
        y = it->second.y();
        r = std::sqrt(x * x + y * y);
        if (r < m_FCal_rmin[imap - 1])
          m_FCal_rmin[imap - 1] = r;
        if (r > m_FCal_rmax[imap - 1])
          m_FCal_rmax[imap - 1] = r;
      }
  }

  auto getClosestFCalCellIndex(
      int layer, float x, float y, int& ieta, int& iphi, int* steps) const
      -> bool
  {
    double rmin = m_FCal_rmin[layer - 21];
    double rmax = m_FCal_rmax[layer - 21];
    int isam = layer - 20;
    double a = 1.;
    const double b = 0.01;
    const int nmax = 100;
    int i = 0;

    const double r = sqrt(x * x + y * y);
    if (r == 0.)
      return false;
    const double r_inverse = 1. / r;

    if ((r / rmax) > (rmin * r_inverse)) {
      x = x * rmax * r_inverse;
      y = y * rmax * r_inverse;
      while ((!m_FCal_ChannelMap.getTileID(isam, a * x, a * y, ieta, iphi))
             && i < nmax)
      {
        a -= b;
        i++;
      }
    } else {
      x = x * rmin * r_inverse;
      y = y * rmin * r_inverse;
      while ((!m_FCal_ChannelMap.getTileID(isam, a * x, a * y, ieta, iphi))
             && i < nmax)
      {
        a += b;
        i++;
      }
    }
    if (steps)
      *steps = i + 1;
    return i < nmax;
  }

  /// @brief Returns the best matching FCAL cell for a hit in a specific layer
  auto get_fcal_cell_id(int layer,
                        double x,
                        double y,
                        double z,
                        int* steps = 0) const -> long long
  {
    int isam = layer - 20;
    int iphi(-100000), ieta(-100000);
    Long64_t mask1[] {0x34, 0x34, 0x35};
    Long64_t mask2[] {0x36, 0x36, 0x37};
    bool found = m_FCal_ChannelMap.getTileID(isam, x, y, ieta, iphi);
    if (steps && found)
      *steps = 0;
    if (!found) {
      // Hit is not matched with any FCal cell! Looking for
      // the closest cell
      found = getClosestFCalCellIndex(layer, x, y, ieta, iphi, steps);
    }
    if (!found) {
      std::cout << "Error: Unable to find the closest FCal cell!" << std::endl;
      return -1;
    }
    Long64_t id = (ieta << 5) + 2 * iphi;
    if (isam == 2)
      id += (8 << 8);

    if (z > 0)
      id += (mask2[isam - 1] << 12);
    else
      id += (mask1[isam - 1] << 12);

    id = id << 44;

    return id;
  }

protected:
  // The constructor argument is set according to
  // https://atlas-sw-doxygen.web.cern.ch/atlas-sw-doxygen/atlas_main--Doxygen/docs/html/d6/d72/FCALChannelMapBuilder_8cxx_source.html
  FCAL_ChannelMap m_FCal_ChannelMap {0};
  std::vector<double> m_FCal_rmin, m_FCal_rmax;
};
