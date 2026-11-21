// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#pragma once

#include <cstddef>
#include <cstdint>
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
  // Base class of alternative geometry handlers (e.g. the FCal)
  virtual ~CaloGeo() = default;
  CaloGeo(const CaloGeo&) = delete;
  auto operator=(const CaloGeo&) -> CaloGeo& = delete;
  CaloGeo(CaloGeo&&) = default;
  auto operator=(CaloGeo&&) -> CaloGeo& = default;

  // Method to build geometry: writes the cell store into the directory
  void build(ROOT::RDataFrame& geo, const std::string& base_path);

  // Method to load geometry from the directory written by build()
  // (rtree_cache_size is unused and only kept for compatibility: the R-trees
  // are built in memory)
  void load(const std::string& base_path,
            size_t rtree_cache_size = std::size_t {5} * 1024 * 1024);

  // Retrieve the best matching cell for a given position
  // Alternative geo handlers need to override this method
  virtual auto get_cell(std::size_t layer, const Position& pos) const
      -> const Cell&;

  // Retrieve the id of the best matching cell for a given position
  auto get_cell_id(std::size_t layer, const Position& pos) const
      -> unsigned long long;
  /// @brief Get the number of cells in a layer
  auto n_cells(std::size_t layer) const -> std::size_t;
  /// @brief Get the total number of cells
  auto n_cells() const -> std::size_t;
  /// @brief Get the number of layers
  auto n_layers() const -> std::size_t;

  /// @brief Get a cell by its ID
  auto get_cell(unsigned long long id) const -> const Cell&;
  /// @brief Get a cell by its index in a layer
  auto get_cell_at_idx(std::size_t layer, size_t idx) const -> const Cell&;

  /// @brief Check if a layer is a barrel layer
  auto is_barrel(std::size_t layer) const -> bool;
  /// @brief Check if cuboid cells in a layer are described by (x, y, z)
  auto is_xyz(std::size_t layer) const -> bool;
  /// @brief Check if cuboid cells in a layer are described by (eta, phi, r)
  auto is_eta_phi_r(std::size_t layer) const -> bool;
  /// @brief Check if cuboid cells in a layer are described by (eta, phi, z)
  auto is_eta_phi_z(std::size_t layer) const -> bool;

  /// @brief Get the z position of a cell in a layer at a given position
  auto zpos(std::size_t layer,
            const Position& pos,
            FastCaloSim::SubPos subpos) const -> double;
  /// @brief Get the r position of a cell in a layer at a given position
  auto rpos(std::size_t layer,
            const Position& pos,
            FastCaloSim::SubPos subpos) const -> double;

  /// @brief Get the min eta extent of a layer on a given detector side
  auto min_eta(std::size_t layer, DetectorSide side) const -> double;
  /// @brief Get the max eta extent of a layer on a given detector side
  auto max_eta(std::size_t layer, DetectorSide side) const -> double;

  /// @brief Set an alternative geometry handler for a layer
  void set_alt_geo_handler(std::size_t layer, std::shared_ptr<CaloGeo> handler);

  /// @brief Set an alternative geometry handler for a range of layers
  void set_alt_geo_handler(std::size_t first_layer,
                           std::size_t last_layer,
                           const std::shared_ptr<CaloGeo>& handler);

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
      eta_extensions[DetectorSide::kEtaPositive] =
          EtaExtremes {.eta_min = 1000, .eta_max = -1000};
      eta_extensions[DetectorSide::kEtaNegative] =
          EtaExtremes {.eta_min = 1000, .eta_max = -1000};
    }
  };

  /// @brief The number of layers
  std::size_t m_n_layers {};
  /// @brief The total number of cells
  std::size_t m_n_total_cells {};

  /// @brief Maps layer id -> R-tree for the cell lookup
  std::unordered_map<std::size_t, std::unique_ptr<RTree>> m_layer_rtrees;

  /// @brief The memory-mapped cell store
  CellStore m_cell_store;

  /// @brief Maps layer ID to vector of cell IDs in that layer
  std::unordered_map<std::size_t, std::vector<unsigned long long>>
      m_layer_cell_ids;

  /// @brief Alternative geo handlers
  /// Allows to implement custom geo handling for specific layers
  /// This is especially useful if the geometry is not well described
  /// by cuboid cells and you want to implement your own cell lookup
  std::unordered_map<std::size_t, std::shared_ptr<CaloGeo>> m_alt_geo_handlers;

  /// @brief Maps layer id -> layer properties
  std::map<std::size_t, LayerFlags> m_layer_flags;

  // What get_cell(layer, pos) needs for a layer, in a single map, so that the
  // lookup done for every hit needs one hash lookup for any layer ID
  struct LayerLookup
  {
    // Alternative geometry handler of the layer (owned by m_alt_geo_handlers)
    const CaloGeo* alt_geo_handler {nullptr};
    // R-tree of the layer (owned by m_layer_rtrees)
    const RTree* rtree {nullptr};
    // Cell store index of each cell of the layer, in R-tree order
    std::vector<std::uint32_t> store_index;
  };
  std::unordered_map<std::size_t, LayerLookup> m_layer_lookup;

  /// @brief Record a cell in the geometry
  void record_cell(std::unique_ptr<Cell> cell);

  /// @brief Update the eta extremes of a layer
  void update_eta_extremes(std::size_t layer, const Cell& cell);
};
