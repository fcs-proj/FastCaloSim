#pragma once

#include <iostream>
#include <map>
#include <vector>

#include <ROOT/RDataFrame.hxx>

#include "FastCaloSim/Geometry/Cell.h"
#include "FastCaloSim/Geometry/FSmap.h"
#include "FastCaloSim/Geometry/RTree.h"

class CaloGeo
{
  using RzMapDouble = std::vector<std::vector<FSmap<double, double>>>;
  using RzMapInt = std::vector<std::vector<FSmap<double, int>>>;

public:
  enum CaloSubPos
  {
    SUBPOS_MID = 0,  // middle
    SUBPOS_ENT = 1,  // entrance
    SUBPOS_EXT = 2  // exit
  };
  explicit CaloGeo(ROOT::RDataFrame& geo) { build(geo); }

  /// @brief Returns the number of cells in a layer
  auto n_cells(int layer) const { return m_layer_tree_map.at(layer).size(); }

  /// @brief Returns the total number of cells in the detector
  auto n_cells() const { return m_n_cells; }

  /// @brief Returns the total number of layers in the detector
  auto n_layers() const { return m_is_barrel_map.size(); }

  /// @brief Returns the cell with a specific cell ID
  auto get_cell(long long id) const -> Cell { return m_cell_id_map.at(id); }

  /// @brief Returns the cell at a specific index in a layer
  auto get_cell_at_idx(int layer, size_t idx) -> Cell
  {
    return m_layer_tree_map.at(layer).at(idx);
  }

  // @brief Returns the best matching cell for a hit in a specific layer
  template<typename T>
  auto get_cell(int layer, const T& hit) const -> Cell
  {
    return m_layer_tree_map.at(layer).query_point(hit);
  }

  // @brief Returns a boolean indicating if a layer is a barrel or endcap
  auto is_barrel(int layer) const { return m_is_barrel_map.at(layer); }

  /// @brief Radial position of the entrance, middle, or exit of a layer at a
  /// specific eta value
  auto rpos(int layer, double eta, CaloSubPos pos) const
  {
    int side = (eta > 0) ? 1 : 0;

    if (pos == CaloSubPos::SUBPOS_ENT) {
      return m_rent_map[side][layer].find_closest(eta)->second;
    }
    if (pos == CaloSubPos::SUBPOS_EXT) {
      return m_rext_map[side][layer].find_closest(eta)->second;
    }

    return m_rmid_map[side][layer].find_closest(eta)->second;
  }

  /// @brief Z position of the entrance, middle, or exit of a layer at a
  /// specific eta value
  auto zpos(int layer, double eta, CaloSubPos pos) const
  {
    int side = (eta > 0) ? 1 : 0;

    if (pos == CaloSubPos::SUBPOS_ENT) {
      return m_zent_map[side][layer].find_closest(eta)->second;
    }
    if (pos == CaloSubPos::SUBPOS_EXT) {
      return m_zext_map[side][layer].find_closest(eta)->second;
    }

    return m_zmid_map[side][layer].find_closest(eta)->second;
  }

  /// @brief Minimum eta extension of the layer
  /// @param layer The layer number
  /// @param eta The eta value that just decides the side of the layer
  auto min_eta(int layer, double eta) const -> double
  {
    int side = (eta > 0) ? 1 : 0;
    return m_min_eta_layer[side][layer];
  }

  /// @brief Maximum eta extension of the layer
  /// @param layer The layer number
  /// @param eta The eta value that just forces the side of the layer
  auto max_eta(int layer, double eta) const -> double
  {
    int side = (eta > 0) ? 1 : 0;
    return m_max_eta_layer[side][layer];
  }

private:
  int m_n_cells {};
  std::map<int, RTree> m_layer_tree_map;
  std::map<long long, Cell> m_cell_id_map;
  std::map<int, bool> m_is_barrel_map;

  std::vector<std::vector<double>> m_min_eta_layer, m_max_eta_layer;
  std::vector<FSmap<double, double>> m_rent_map[2], m_zent_map[2],
      m_rmid_map[2], m_zmid_map[2], m_rext_map[2], m_zext_map[2];

  void record_cell(const Cell& cell)
  {
    m_layer_tree_map[cell->layer()].insert_cell(cell);
    m_cell_id_map.emplace(cell->id(), cell);
  }

  void build(ROOT::RDataFrame& geo)
  {
    m_n_cells = *geo.Count();

    auto layer = geo.Take<long long>("layer");
    auto isBarrel = geo.Take<long long>("isBarrel");
    auto id = geo.Take<long long>("id");
    auto x = geo.Take<double>("x");
    auto y = geo.Take<double>("y");
    auto z = geo.Take<double>("z");
    auto phi = geo.Take<double>("phi");
    auto eta = geo.Take<double>("eta");
    auto r = geo.Take<double>("r");
    auto dx = geo.Take<double>("dx");
    auto dy = geo.Take<double>("dy");
    auto dz = geo.Take<double>("dz");
    auto dphi = geo.Take<double>("dphi");
    auto deta = geo.Take<double>("deta");
    auto dr = geo.Take<double>("dr");
    auto isCartesian = geo.Take<long long>("isCartesian");
    auto isCylindrical = geo.Take<long long>("isCylindrical");
    auto isECCylindrical = geo.Take<long long>("isECCylindrical");

    for (size_t i = 0; i < m_n_cells; ++i) {
      m_is_barrel_map.emplace(layer->at(i), static_cast<bool>(isBarrel->at(i)));

      CaloPos caloPos {
          x->at(i), y->at(i), z->at(i), eta->at(i), phi->at(i), r->at(i)};

      if (isCartesian->at(i)) {
        record_cell(XYZCell {id->at(i),
                             caloPos,
                             layer->at(i),
                             static_cast<bool>(isBarrel->at(i)),
                             dx->at(i),
                             dy->at(i),
                             dz->at(i)});
      } else if (isCylindrical->at(i)) {
        record_cell(EtaPhiRCell {id->at(i),
                                 caloPos,
                                 layer->at(i),
                                 static_cast<bool>(isBarrel->at(i)),
                                 deta->at(i),
                                 dphi->at(i),
                                 dr->at(i)});
      } else if (isECCylindrical->at(i)) {
        record_cell(EtaPhiZCell {id->at(i),
                                 caloPos,
                                 layer->at(i),
                                 static_cast<bool>(isBarrel->at(i)),
                                 deta->at(i),
                                 dphi->at(i),
                                 dz->at(i)});
      }
    }

    init_rz_maps();
  }

  void init_rz_maps()
  {
    /**
     * @brief Initializes the RZ maps for the detector.
     *
     * This method initializes and populates several maps used to aggregate and
     * average various properties of detector cells, organized by detector side
     * and layer. These properties include the pseudorapidity (eta), radial
     * position (r), and axial position (z), and they essentially encode the
     * entrance, mid, and exit positions of the calorimeter layers.
     *
     * The method performs the following steps:
     * 1. Initializes maps to store accumulated sums and counts of cell
     * properties:
     * 2. Iterates over each cell in each layer to accumulate these properties
     *    in the corresponding maps.
     * 3. Computes averages for each property by dividing the accumulated sums
     *    by the counts and stores the results in member maps for further
     * analysis.
     */

    int MAX_LAYER = n_layers();
    RzMapDouble rz_map_eta(2, std::vector<FSmap<double, double>>(MAX_LAYER));
    RzMapDouble rz_map_rmid(2, std::vector<FSmap<double, double>>(MAX_LAYER));
    RzMapDouble rz_map_zmid(2, std::vector<FSmap<double, double>>(MAX_LAYER));
    RzMapDouble rz_map_rent(2, std::vector<FSmap<double, double>>(MAX_LAYER));
    RzMapDouble rz_map_zent(2, std::vector<FSmap<double, double>>(MAX_LAYER));
    RzMapDouble rz_map_rext(2, std::vector<FSmap<double, double>>(MAX_LAYER));
    RzMapDouble rz_map_zext(2, std::vector<FSmap<double, double>>(MAX_LAYER));
    RzMapInt rz_map_n(2, std::vector<FSmap<double, int>>(MAX_LAYER));

    // Resize for two sides
    m_min_eta_layer.resize(2);
    m_max_eta_layer.resize(2);

    for (int side = 0; side <= 1; ++side) {
      // Resize for MAX_LAYER layers
      m_min_eta_layer[side].resize(MAX_LAYER, +1000);
      m_max_eta_layer[side].resize(MAX_LAYER, -1000);

      m_rent_map[side].resize(MAX_LAYER);
      m_zent_map[side].resize(MAX_LAYER);
      m_rmid_map[side].resize(MAX_LAYER);
      m_zmid_map[side].resize(MAX_LAYER);
      m_rext_map[side].resize(MAX_LAYER);
      m_zext_map[side].resize(MAX_LAYER);

      // Initialize min and max eta values for each layer
      for (int layer = 0; layer < MAX_LAYER; ++layer) {
        m_min_eta_layer.at(side).at(layer) = +1000;
        m_max_eta_layer.at(side).at(layer) = -1000;
      }
    }

    for (const auto& [layer, cells] : m_layer_tree_map) {
      for (size_t i = 0; i < cells.size(); ++i) {
        const auto& cell = cells.at(i);
        int side = (cell->eta() > 0) ? 1 : 0;
        int sign_side = (side == 1) ? +1 : -1;
        double eta = cell->eta();

        double deta = cell->isXYZ() ? 0 : cell->deta();

        double min_eta = cell->eta() - deta / 2;
        double max_eta = cell->eta() + deta / 2;

        m_min_eta_layer[side][layer] =
            std::min(m_min_eta_layer[side][layer], min_eta);
        m_max_eta_layer[side][layer] =
            std::max(m_max_eta_layer[side][layer], max_eta);

        rz_map_eta[side][layer][eta] += cell->eta();
        rz_map_rmid[side][layer][eta] += cell->r();
        rz_map_zmid[side][layer][eta] += cell->z();

        double dz = cell->isEtaPhiR() ? 0 : cell->dz();
        double dr = cell->isEtaPhiR() ? cell->dr() : 0;

        rz_map_zent[side][layer][eta] += cell->z() - dz / 2 * sign_side;
        rz_map_zext[side][layer][eta] += cell->z() + dz / 2 * sign_side;

        rz_map_rent[side][layer][eta] += cell->r() - dr / 2;
        rz_map_rext[side][layer][eta] += cell->r() + dr / 2;

        rz_map_n[side][layer][eta]++;
      }
    }

    for (int side = 0; side <= 1; ++side) {
      for (int layer = 0; layer < MAX_LAYER; ++layer) {
        if (!rz_map_n[side][layer].empty()) {
          for (const auto& [eta, count] : rz_map_n[side][layer]) {
            double avg_eta = rz_map_eta[side][layer][eta] / count;
            double rmid = rz_map_rmid[side][layer][eta] / count;
            double zmid = rz_map_zmid[side][layer][eta] / count;
            double rent = rz_map_rent[side][layer][eta] / count;
            double zent = rz_map_zent[side][layer][eta] / count;
            double rext = rz_map_rext[side][layer][eta] / count;
            double zext = rz_map_zext[side][layer][eta] / count;

            m_rent_map[side][layer][avg_eta] = rent;
            m_zent_map[side][layer][avg_eta] = zent;
            m_rmid_map[side][layer][avg_eta] = rmid;
            m_zmid_map[side][layer][avg_eta] = zmid;
            m_rext_map[side][layer][avg_eta] = rext;
            m_zext_map[side][layer][avg_eta] = zext;
          }
        } else {
          std::cout << "rz-map for side=" << side << " layer=" << layer
                    << " is empty!!!" << std::endl;
        }
      }
    }
  }
};
