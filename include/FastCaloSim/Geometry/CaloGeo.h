// Copyright (c) 2024 CERN for the benefit of the FastCaloSim project

#pragma once

#include <map>

#include <FastCaloSim/FastCaloSim_export.h>

#include <ROOT/RDataFrame.hxx>

#include "FastCaloSim/Geometry/Cell.h"
#include "FastCaloSim/Geometry/RTree.h"

class FASTCALOSIM_EXPORT CaloGeo
{
public:
  enum DetectorSide
  {
    kEtaPositive,
    kEtaNegative
  };

  // Default constructor
  CaloGeo() = default;

  // Constructor from a ROOT dataframe
  explicit CaloGeo(ROOT::RDataFrame& geo);

  // Retrieve the id of the best matching cell for a given position
  // Alternative geo handlers need to override this method
  virtual auto get_cell_id(int layer, const Position& pos) const -> long long;

  // Retrieve the best matching cell for a given position
  auto get_cell(int layer, const Position& pos) const -> Cell;

  /// @brief Get the number of cells in a layer
  auto n_cells(int layer) const -> unsigned int;
  /// @brief Get the total number of cells
  auto n_cells() const -> unsigned int;
  /// @brief Get the number of layers
  auto n_layers() const -> unsigned int;

  /// @brief Get a cell by its ID
  auto get_cell(long long id) const -> Cell;
  /// @brief Get a cell by its index in a layer
  auto get_cell_at_idx(int layer, size_t idx) -> Cell;

  /// @brief Check if a layer is a barrel layer
  auto is_barrel(int layer) const -> bool;
  /// @brief Check if cuboid cells in a layer are described by (x, y, z)
  auto is_xyz(int layer) const -> bool;
  /// @brief Check if cuboid cells in a layer are described by (eta, phi, r)
  auto is_eta_phi_r(int layer) const -> bool;
  /// @brief Check if cuboid cells in a layer are described by (eta, phi, z)
  auto is_eta_phi_z(int layer) const -> bool;

  /// @brief Get the z position of a cell in a layer at a given position
  auto zpos(int layer,
            const Position& pos,
            Cell::SubPos subpos) const -> double;
  /// @brief Get the r position of a cell in a layer at a given position
  auto rpos(int layer,
            const Position& pos,
            Cell::SubPos subpos) const -> double;

  /// @brief Get the min eta extent of a layer on a given detector side
  auto min_eta(int layer, DetectorSide side) const -> double;
  /// @brief Get the max eta extent of a layer on a given detector side
  auto max_eta(int layer, DetectorSide side) const -> double;

  /// @brief Set an alternative geometry handler for a layer
  void set_alt_geo_handler(int layer, std::shared_ptr<CaloGeo> handler);

  /// @brief Set an alternative geometry handler for a range of layers
  void set_alt_geo_handler(int first_layer,
                           int last_layer,
                           std::shared_ptr<CaloGeo> handler);

private:
  /// @brief The eta extent of a layer on a given detector side
  struct EtaExtremes
  {
    double eta_min, eta_max;
  };

  /// @brief Properties of a layer
  struct LayerFlags
  {
    // Is the layer a barrel layer?
    bool is_barrel {false};
    // How are the cells in the layer described?
    bool is_xyz {false};
    bool is_eta_phi_r {false};
    bool is_eta_phi_z {false};
    bool is_r_phi_z {false};
    // For each detector side, what are the eta_min and eta_max values?
    std::map<DetectorSide, EtaExtremes> eta_extensions;

    // Define default values for the eta extensions of the layers
    LayerFlags()
    {
      eta_extensions[DetectorSide::kEtaPositive] = EtaExtremes {1000, -1000};
      eta_extensions[DetectorSide::kEtaNegative] = EtaExtremes {1000, -1000};
    }
  };

  /// @brief The number of layers
  int m_n_layers {};
  /// @brief The total number of cells
  int m_n_total_cells {};
  /// @brief Maps layer id -> RTree
  std::unordered_map<int, RTree> m_layer_tree_map;
  /// @brief Maps cell id -> Cell
  std::unordered_map<long long, Cell> m_cell_id_map;

  /// @brief Alternative geo handlers
  /// Allows to implement custom geo handling for specific layers
  /// This is especially useful if the geometry is not well described
  /// by cuboid cells and you want to implement your own cell lookup
  std::unordered_map<int, std::shared_ptr<CaloGeo>> m_alt_geo_handlers;

  /// @brief Maps layer id -> layer properties
  std::map<int, LayerFlags> m_layer_flags;
  /// @brief Record a cell in the geometry
  void record_cell(const Cell& cell);
  /// @brief Update the eta extremes of a layer
  void update_eta_extremes(int layer, const Cell& cell);
  /// @brief Build the geometry from a ROOT dataframe
  void build(ROOT::RDataFrame& geo);
};
