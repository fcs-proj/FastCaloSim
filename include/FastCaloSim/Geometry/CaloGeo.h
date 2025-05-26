// Copyright (c) 2024 CERN for the benefit of the FastCaloSim project

#pragma once

#include <map>
#include <memory>
#include <unordered_map>
#include <vector>

#include <FastCaloSim/FastCaloSim_export.h>

#include <ROOT/RDataFrame.hxx>

#include "FastCaloSim/Geometry/Cell.h"
#include "FastCaloSim/Geometry/CellStore.h"
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

  // Method to build geometry
  void build(ROOT::RDataFrame& geo, const std::string& rtree_base_path);

  // Method to load geometry
  void load(const std::string& rtree_base_path,
            size_t rtree_cache_size = 5 * 1024 * 1024,
            size_t cell_store_cache_size = 10 * 1024 * 1024);

  // Retrieve the id of the best matching cell for a given position
  // Alternative geo handlers need to override this method
  virtual auto get_cell_id(unsigned int layer,
                           const Position& pos) const -> unsigned long long;

  // Retrieve the best matching cell for a given position
  virtual auto get_cell(unsigned int layer,
                        const Position& pos) const -> const Cell&;

  /// @brief Get the number of cells in a layer
  auto n_cells(unsigned int layer) const -> unsigned int;
  /// @brief Get the total number of cells
  auto n_cells() const -> unsigned int;
  /// @brief Get the number of layers
  auto n_layers() const -> unsigned int;

  /// @brief Get a cell by its ID
  auto get_cell(unsigned long long id) const -> const Cell&;
  /// @brief Get a cell by its index in a layer
  auto get_cell_at_idx(unsigned int layer, size_t idx) const -> const Cell&;

  /// @brief Check if a layer is a barrel layer
  auto is_barrel(unsigned int layer) const -> bool;
  /// @brief Check if cuboid cells in a layer are described by (x, y, z)
  auto is_xyz(unsigned int layer) const -> bool;
  /// @brief Check if cuboid cells in a layer are described by (eta, phi, r)
  auto is_eta_phi_r(unsigned int layer) const -> bool;
  /// @brief Check if cuboid cells in a layer are described by (eta, phi, z)
  auto is_eta_phi_z(unsigned int layer) const -> bool;

  /// @brief Get the z position of a cell in a layer at a given position
  auto zpos(unsigned int layer,
            const Position& pos,
            Cell::SubPos subpos) const -> double;
  /// @brief Get the r position of a cell in a layer at a given position
  auto rpos(unsigned int layer,
            const Position& pos,
            Cell::SubPos subpos) const -> double;

  /// @brief Get the min eta extent of a layer on a given detector side
  auto min_eta(unsigned int layer, DetectorSide side) const -> double;
  /// @brief Get the max eta extent of a layer on a given detector side
  auto max_eta(unsigned int layer, DetectorSide side) const -> double;

  /// @brief Set an alternative geometry handler for a layer
  void set_alt_geo_handler(unsigned int layer,
                           std::shared_ptr<CaloGeo> handler);

  /// @brief Set an alternative geometry handler for a range of layers
  void set_alt_geo_handler(unsigned int first_layer,
                           unsigned int last_layer,
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
    RTreeHelpers::CoordinateSystem coordinate_system {
        RTreeHelpers::CoordinateSystem::Undefined};
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
  unsigned int m_n_layers {};
  /// @brief The total number of cells
  unsigned int m_n_total_cells {};

  /// @brief Maps layer id -> RTreeQuery (for loaded trees)
  std::unordered_map<unsigned int, std::unique_ptr<RTreeQuery>>
      m_layer_rtree_queries;

  /// @brief The memory-mapped cell store
  CellStore m_cell_store;

  /// @brief Maps layer ID to vector of cell IDs in that layer
  std::unordered_map<unsigned int, std::vector<unsigned long long>>
      m_layer_cell_ids;

  /// @brief Alternative geo handlers
  /// Allows to implement custom geo handling for specific layers
  /// This is especially useful if the geometry is not well described
  /// by cuboid cells and you want to implement your own cell lookup
  std::unordered_map<unsigned int, std::shared_ptr<CaloGeo>> m_alt_geo_handlers;

  /// @brief Maps layer id -> layer properties
  std::map<unsigned int, LayerFlags> m_layer_flags;

  /// @brief Record a cell in the geometry
  void record_cell(std::unique_ptr<Cell> cell);

  /// @brief Update the eta extremes of a layer
  void update_eta_extremes(unsigned int layer, const Cell& cell);
};
