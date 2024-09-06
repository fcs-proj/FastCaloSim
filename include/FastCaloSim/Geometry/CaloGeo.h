#pragma once

#include <map>
#include <vector>

#include <ROOT/RDataFrame.hxx>

#include "FastCaloSim/Geometry/Cell.h"
#include "FastCaloSim/Geometry/RTree.h"

class CaloGeo
{
public:
  /// @brief Enum to indicate the side of the detector η > 0 or η < 0
  enum DetectorSide
  {
    kEtaPositive,
    kEtaNegative
  };

  explicit CaloGeo(ROOT::RDataFrame& geo) { build(geo); }

  /// @brief Returns the number of cells in a layer
  auto n_cells(int layer) const { return m_layer_tree_map.at(layer).size(); }

  /// @brief Returns the total number of cells in the detector
  auto n_cells() const { return m_n_total_cells; }

  /// @brief Returns the total number of layers in the detector
  auto n_layers() const { return m_n_layers; }

  /// @brief Returns the cell with a specific cell ID
  auto get_cell(long long id) const -> Cell { return m_cell_id_map.at(id); }

  /// @brief Returns the cell at a specific index in a layer
  auto get_cell_at_idx(int layer, size_t idx) -> Cell
  {
    return m_layer_tree_map.at(layer).at(idx);
  }

  /// @brief Returns the best matching cell for a hit in a specific layer
  // Note that the hit must have the same coordinate system as the cells in the
  // tree
  template<typename T>
  auto get_cell(int layer, const T& hit) const -> Cell
  {
    return m_layer_tree_map.at(layer).query_point(hit);
  }

  /// @brief Returns a boolean indicating if a layer is a barrel or endcap
  auto is_barrel(int layer) const { return m_layer_flags.at(layer).is_barrel; };

  /// @brief Returns a boolean indicating if a layer is in XYZ coordinates
  auto is_xyz(int layer) const { return m_layer_flags.at(layer).ix_xyz; };

  /// @brief Returns a boolean indicating if a layer is in EtaPhiR coordinates
  auto is_eta_phi_r(int layer) const
  {
    return m_layer_flags.at(layer).is_eta_phi_r;
  };

  /// @brief Returns a boolean indicating if a layer is in EtaPhiZ coordinates
  auto is_eta_phi_z(int layer) const
  {
    return m_layer_flags.at(layer).is_eta_phi_z;
  };

  /// @brief Z position of the entrance, middle, or exit of a layer for
  // the cell closest to a hit (either in x,y or in eta, phi space)
  template<typename T>
  auto zpos(int layer, const T& hit, Cell::SubPos subpos) const -> double
  {
    const auto& cell = get_cell(layer, hit);

    return cell.z(subpos);
  }

  /// @brief R position of the entrance, middle, or exit of a layer for
  /// the cell closest to a hit (either in x,y or in eta, phi space)
  template<typename T>
  auto rpos(int layer, const T& hit, Cell::SubPos subpos) const -> double
  {
    const auto& cell = get_cell(layer, hit);

    return cell.r(subpos);
  }

  /// @brief Returns the minmax eta extension of a layer for a side of the
  // detector
  auto min_eta(int layer, DetectorSide side) const
  {
    return m_layer_flags.at(layer).eta_extensions.at(side).eta_min;
  }

  /// @brief Returns the max eta extension of a layer for a side of the detector
  auto max_eta(int layer, DetectorSide side) const
  {
    return m_layer_flags.at(layer).eta_extensions.at(side).eta_max;
  }

private:
  /// @brief Struct to store the minimum and maximum eta extensions of a layer
  struct EtaExtremes
  {
    double eta_min, eta_max;
  };

  /// @brief Struct to store properties  of a layer
  struct LayerFlags
  {
    bool is_barrel {false};
    bool ix_xyz {false};
    bool is_eta_phi_r {false};
    bool is_eta_phi_z {false};
    /// @brief The eta extensions for the positive and negative detector half
    // spaces
    std::map<DetectorSide, EtaExtremes> eta_extensions;

    // Define default values for the eta extensions of the layers
    LayerFlags()
    {
      eta_extensions[DetectorSide::kEtaPositive] = EtaExtremes {1000, -1000};
      eta_extensions[DetectorSide::kEtaNegative] = EtaExtremes {1000, -1000};
    }
  };

  int m_n_layers {};
  int m_n_total_cells {};
  std::map<int, RTree> m_layer_tree_map;
  std::map<long long, Cell> m_cell_id_map;
  std::unordered_map<int, LayerFlags> m_layer_flags;

  /// @brief Record the cell in the geometry
  void record_cell(const Cell& cell)
  {
    // Insert cell RTree
    auto layer = cell.layer();
    m_layer_tree_map[layer].insert_cell(cell);
    // Insert cell in cell ID -> cell map
    m_cell_id_map.emplace(cell.id(), cell);
  }

  /// @brief For each layer, we pre-compute the maximum
  // and minimum eta extensions of the layers for
  // the positive (η > 0) and negative (η < 0) detector half spaces
  void update_eta_extremes(int layer, const Cell& cell)
  {
    // Choose detector side
    DetectorSide side = cell.eta() > 0 ? kEtaPositive : kEtaNegative;

    // Cell half width
    // for FCAL we don't take into account the half width
    double half_width = cell.isXYZ() ? 0 : cell.deta() / 2;

    // Compute eta extremes for the side
    double min_eta = cell.eta() - half_width;
    double max_eta = cell.eta() + half_width;

    auto& eta_extremes = m_layer_flags.at(layer).eta_extensions.at(side);

    // Update the eta extremes
    if (max_eta > eta_extremes.eta_max) {
      eta_extremes.eta_max = max_eta;
    }
    if (min_eta < eta_extremes.eta_min) {
      eta_extremes.eta_min = min_eta;
    }
  }

  /// @brief Build the geometry from the RDataFrame
  void build(ROOT::RDataFrame& geo)
  {
    // Record number of layers
    m_n_layers = geo.Max<long long>("layer").GetValue();

    // Initialize maximum and minimum eta extensions of the layers
    for (int i = 0; i < m_n_layers + 1; ++i) {
      m_layer_flags.emplace(i, LayerFlags {});
    }

    // Record number of total cells
    m_n_total_cells = *geo.Count();

    // Record the cell information
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

    // TODO: Currently layer, isBarrel, isCartesian, isCylindrical, and
    // isECCylindrical are recorded by ROOT as long long and require casting We
    // should likely change how the types are written in the future

    // Loop over all cells in the detector
    for (size_t i = 0; i < m_n_total_cells; ++i) {
      bool is_barrel = static_cast<bool>(isBarrel->at(i));
      bool is_XYZ = static_cast<bool>(isCartesian->at(i));
      bool is_EtaPhiR = static_cast<bool>(isCylindrical->at(i));
      bool is_EtaPhiZ = static_cast<bool>(isECCylindrical->at(i));

      m_layer_flags.at(layer->at(i)).is_barrel = is_barrel;
      m_layer_flags.at(layer->at(i)).ix_xyz = is_XYZ;
      m_layer_flags.at(layer->at(i)).is_eta_phi_r = is_EtaPhiR;
      m_layer_flags.at(layer->at(i)).is_eta_phi_z = is_EtaPhiZ;

      // Position of the center of the cell
      Position pos = Position {.m_x = x->at(i),
                               .m_y = y->at(i),
                               .m_z = z->at(i),
                               .m_eta = eta->at(i),
                               .m_phi = phi->at(i),
                               .m_r = r->at(i)};
      // Construct the cell
      Cell cell = Cell {id->at(i),
                        pos,
                        layer->at(i),
                        is_barrel,
                        is_XYZ,
                        is_EtaPhiR,
                        is_EtaPhiZ,
                        dx->at(i),
                        dy->at(i),
                        dz->at(i),
                        deta->at(i),
                        dphi->at(i),
                        dr->at(i)};

      // Record the cell in the geometry
      record_cell(cell);
      // Update the edge cells
      update_eta_extremes(layer->at(i), cell);
    }
  }
};
