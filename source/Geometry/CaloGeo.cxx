// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <format>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>

#include "FastCaloSim/Geometry/CaloGeo.h"

#include "FastCaloSim/Core/MLogging.h"
#include "FastCaloSim/Geometry/CellStoreBuilder.h"

namespace
{
// The cell lookup and the hit-cell distance work in (eta, phi) or (x, y);
// cells described in (r, phi, z) have no eta extent. Checked for every cell,
// not only the first of a layer.
void reject_r_phi_z_cell(bool is_r_phi_z, std::size_t layer_id)
{
  if (is_r_phi_z) {
    throw std::runtime_error(std::format(
        "Layer {}: cells in (r, phi, z) coordinates (isRPhiZ) are not "
        "supported by the cell lookup",
        layer_id));
  }
}

// How the cells of a layer are described, from the flags of one of its cells
auto coordinate_system_of(bool is_xyz, bool is_eta_phi_r, bool is_eta_phi_z)
    -> RTreeHelpers::CoordinateSystem
{
  if (is_xyz) {
    return RTreeHelpers::CoordinateSystem::XYZ;
  }
  if (is_eta_phi_r) {
    return RTreeHelpers::CoordinateSystem::EtaPhiR;
  }
  if (is_eta_phi_z) {
    return RTreeHelpers::CoordinateSystem::EtaPhiZ;
  }
  return RTreeHelpers::CoordinateSystem::Undefined;
}
// The layer IDs of a geometry, in increasing order
auto sorted_layer_ids(const std::unordered_set<std::size_t>& layer_ids)
    -> std::vector<std::size_t>
{
  std::vector<std::size_t> sorted(layer_ids.begin(), layer_ids.end());
  std::ranges::sort(sorted);
  return sorted;
}
}  // namespace

auto CaloGeo::layer_lookup(std::size_t layer) -> LayerLookup&
{
  if (layer >= n_dense_layer_ids) {
    return m_sparse_layer_lookup[layer];
  }
  if (layer >= m_layer_lookup.size()) {
    m_layer_lookup.resize(layer + 1);
  }
  return m_layer_lookup[layer];
}

auto CaloGeo::get_cell(std::size_t layer, const Position& pos) const -> Cell
{
  if (const LayerLookup* found = find_layer_lookup(layer)) {
    const LayerLookup& lookup = *found;

    // Check if an alternative geometry handler is set for the layer
    if (lookup.alt_geo_handler != nullptr) {
      // Delegate to alternative geometry handler
      Cell cell = lookup.alt_geo_handler->get_cell(layer, pos);
      // Make sure returned cell is valid
      if (!cell.is_valid()) {
        throw std::runtime_error(
            "Invalid cell ID from alternative geometry handler");
      }
      return cell;
    }

    // Default path: R-tree lookup, which gives the cell's index in the store
    if (lookup.rtree != nullptr) {
      return m_cell_store.get_at_index(
          lookup.store_index[lookup.rtree->query_point_index(pos)]);
    }
  }
  throw std::runtime_error("No RTree loaded for layer "
                           + std::to_string(layer));
}

auto CaloGeo::get_cell_id(std::size_t layer, const Position& pos) const
    -> unsigned long long
{
  const Cell cell = get_cell(layer, pos);
  return cell.id();
}

auto CaloGeo::n_cells() const -> std::size_t
{
  return m_n_total_cells;
}

auto CaloGeo::n_cells(std::size_t layer) const -> std::size_t
{
  auto it = m_layer_cell_ids.find(layer);
  if (it != m_layer_cell_ids.end()) {
    return it->second.size();
  }
  return 0;
}

auto CaloGeo::n_layers() const -> std::size_t
{
  return m_layer_ids.size();
}

auto CaloGeo::layer_ids() const -> const std::vector<std::size_t>&
{
  return m_layer_ids;
}

auto CaloGeo::layer_id_bound() const -> std::size_t
{
  return m_layer_ids.empty() ? 0 : m_layer_ids.back() + 1;
}

auto CaloGeo::has_layer(std::size_t layer) const -> bool
{
  return std::ranges::binary_search(m_layer_ids, layer);
}

auto CaloGeo::get_cell(unsigned long long id) const -> Cell
{
  return m_cell_store.get(id);
}

auto CaloGeo::get_cell_at_idx(std::size_t layer, size_t idx) const -> Cell
{
  auto layer_it = m_layer_cell_ids.find(layer);
  if (layer_it == m_layer_cell_ids.end() || idx >= layer_it->second.size()) {
    throw std::runtime_error("Invalid layer ID or index out of bounds");
  }

  unsigned long long const cell_id = layer_it->second[idx];
  return m_cell_store.get(cell_id);
}

auto CaloGeo::is_barrel(std::size_t layer) const -> bool
{
  return m_layer_flags.at(layer).is_barrel;
}

auto CaloGeo::coordinate_system(std::size_t layer) const
    -> RTreeHelpers::CoordinateSystem
{
  const LayerLookup* lookup = find_layer_lookup(layer);
  if (lookup == nullptr || !lookup->coordinate_system) {
    throw std::out_of_range("No cells in layer " + std::to_string(layer));
  }
  return *lookup->coordinate_system;
}

auto CaloGeo::is_xyz(std::size_t layer) const -> bool
{
  return coordinate_system(layer) == RTreeHelpers::CoordinateSystem::XYZ;
}

auto CaloGeo::is_eta_phi_r(std::size_t layer) const -> bool
{
  return coordinate_system(layer) == RTreeHelpers::CoordinateSystem::EtaPhiR;
}

auto CaloGeo::is_eta_phi_z(std::size_t layer) const -> bool
{
  return coordinate_system(layer) == RTreeHelpers::CoordinateSystem::EtaPhiZ;
}

auto CaloGeo::zpos(std::size_t layer,
                   const Position& pos,
                   FastCaloSim::SubPos subpos) const -> double
{
  const Cell cell = get_cell(layer, pos);
  return cell.z(subpos);
}

auto CaloGeo::rpos(std::size_t layer,
                   const Position& pos,
                   FastCaloSim::SubPos subpos) const -> double
{
  const Cell cell = get_cell(layer, pos);
  return cell.r(subpos);
}

auto CaloGeo::min_eta(std::size_t layer, DetectorSide side) const -> double
{
  return m_layer_flags.at(layer).eta_extensions.at(side).eta_min;
}

auto CaloGeo::max_eta(std::size_t layer, DetectorSide side) const -> double
{
  return m_layer_flags.at(layer).eta_extensions.at(side).eta_max;
}

void CaloGeo::set_alt_geo_handler(std::size_t layer,
                                  std::shared_ptr<CaloGeo> handler)
{
  layer_lookup(layer).alt_geo_handler = std::move(handler);
}

void CaloGeo::set_alt_geo_handler(std::size_t first_layer,
                                  std::size_t last_layer,
                                  const std::shared_ptr<CaloGeo>& handler)
{
  for (std::size_t i = first_layer; i <= last_layer; ++i) {
    layer_lookup(i).alt_geo_handler = handler;
  }
}

void CaloGeo::set_layer_roles(FastCaloSim::LayerRoles roles)
{
  m_layer_roles = std::move(roles);
}

auto CaloGeo::has_layer_roles() const -> bool
{
  return m_layer_roles.has_value();
}

void CaloGeo::throw_no_layer_roles()
{
  throw std::runtime_error(
      "CaloGeo: the roles of the layers are not set; the host must set them "
      "with CaloGeo::set_layer_roles()");
}

void CaloGeo::update_eta_extremes(std::size_t layer, const Cell& cell)
{
  DetectorSide const side = cell.eta() > 0 ? kEtaPositive : kEtaNegative;
  double const half_width = cell.deta() / 2;
  double const min_eta = cell.eta() - half_width;
  double const max_eta = cell.eta() + half_width;
  auto& eta_extremes = m_layer_flags.at(layer).eta_extensions.at(side);
  eta_extremes.eta_max = std::max(max_eta, eta_extremes.eta_max);
  eta_extremes.eta_min = std::min(min_eta, eta_extremes.eta_min);
}

// Method to build geometry
void CaloGeo::build(ROOT::RDataFrame& geo, const std::string& base_path)
{
  // Start timing
  auto start_time = std::chrono::high_resolution_clock::now();

  // Pre-load the geometry data into memory
  auto df = geo.Cache(geo.GetColumnNames());

  // Get all column data
  auto layer = df.Take<unsigned int>("layer");
  auto isBarrel = df.Take<bool>("isBarrel");
  auto id = df.Take<unsigned long long>("id");
  auto x = df.Take<float>("x");
  auto y = df.Take<float>("y");
  auto z = df.Take<float>("z");
  auto phi = df.Take<float>("phi");
  auto eta = df.Take<float>("eta");
  auto r = df.Take<float>("r");
  auto dx = df.Take<float>("dx");
  auto dy = df.Take<float>("dy");
  auto dz = df.Take<float>("dz");
  auto dphi = df.Take<float>("dphi");
  auto deta = df.Take<float>("deta");
  auto dr = df.Take<float>("dr");
  auto isXYZ = df.Take<bool>("isXYZ");
  auto isEtaPhiR = df.Take<bool>("isEtaPhiR");
  auto isEtaPhiZ = df.Take<bool>("isEtaPhiZ");
  auto isRPhiZ = df.Take<bool>("isRPhiZ");

  // Count total cells
  m_n_total_cells = layer->size();

  // Find unique layers efficiently
  std::unordered_set<std::size_t> unique_layers;
  for (const auto& layer_id : *layer) {
    unique_layers.insert(layer_id);
  }
  m_layer_ids = sorted_layer_ids(unique_layers);

  // Pre-allocate layer flags with proper initialization for actual layer IDs
  for (const auto& layer_id : unique_layers) {
    m_layer_flags[layer_id] = LayerFlags {};

    // Pre-initialize the eta extensions to extreme values
    m_layer_flags[layer_id].eta_extensions[kEtaPositive].eta_min =
        std::numeric_limits<double>::max();
    m_layer_flags[layer_id].eta_extensions[kEtaPositive].eta_max =
        -std::numeric_limits<double>::max();
    m_layer_flags[layer_id].eta_extensions[kEtaNegative].eta_min =
        std::numeric_limits<double>::max();
    m_layer_flags[layer_id].eta_extensions[kEtaNegative].eta_max =
        -std::numeric_limits<double>::max();
  }

  // Track which layers we've already set flags for to avoid redundant
  std::unordered_set<std::size_t> processed_layer_flags;
  // This will be used to build our memory-mapped cell store
  CellStoreBuilder cell_store_builder;

  // Process cells in bulk
  for (size_t i = 0; i < m_n_total_cells; ++i) {
    unsigned int layer_id = layer->at(i);
    reject_r_phi_z_cell(isRPhiZ->at(i), layer_id);

    // Only update layer flags once per layer
    if (!processed_layer_flags.contains(layer_id)) {
      m_layer_flags[layer_id].is_barrel = isBarrel->at(i);
      layer_lookup(layer_id).coordinate_system = coordinate_system_of(
          isXYZ->at(i), isEtaPhiR->at(i), isEtaPhiZ->at(i));
      processed_layer_flags.insert(layer_id);
    }

    // The geometry is stored in single precision
    const auto value = [i](auto& column) -> auto
    { return static_cast<double>(column->at(i)); };
    const Position pos {.m_x = value(x),
                        .m_y = value(y),
                        .m_z = value(z),
                        .m_eta = value(eta),
                        .m_phi = value(phi),
                        .m_r = value(r)};

    Cell const cell(id->at(i),
                    pos,
                    layer_id,
                    isBarrel->at(i),
                    isXYZ->at(i),
                    isEtaPhiR->at(i),
                    isEtaPhiZ->at(i),
                    isRPhiZ->at(i),
                    value(dx),
                    value(dy),
                    value(dz),
                    value(deta),
                    value(dphi),
                    value(dr));

    update_eta_extremes(layer_id, cell);
    m_layer_cell_ids[layer_id].push_back(cell.id());

    cell_store_builder.add_cell(cell);
  }

  // Write the cell store to disk
  // Memory will be freed after writing
  cell_store_builder.write(base_path + "/cellstore");

  auto end_time = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> const elapsed = end_time - start_time;

  // Print timing information
  // CaloGeo has no output level
  FCS_MSG_LVL_NOCHK(::FastCaloSim::MSG::INFO,
                    "Done building calo geometry. Took "
                        << elapsed.count() << " s (" << m_n_total_cells
                        << " cells in " << n_layers() << " layers)");
}

// Method to load geometry
void CaloGeo::load(const std::string& base_path, size_t /*rtree_cache_size*/)
{
  // Start timing
  auto start_time = std::chrono::high_resolution_clock::now();

  // Invalidate the layer lookups before the cell store is replaced, so that a
  // failed load leaves none that point into the new store
  const auto invalidate = [](LayerLookup& lookup) -> void
  {
    lookup.rtree.reset();
    lookup.store_index.clear();
  };
  std::ranges::for_each(m_layer_lookup, invalidate);
  for (auto& entry : m_sparse_layer_lookup) {
    invalidate(entry.second);
  }

  // Load the cell store
  m_cell_store.load(base_path + "/cellstore");

  // Rebuild the layer cell ID map
  m_layer_cell_ids.clear();
  m_n_total_cells = 0;
  std::unordered_set<std::size_t> unique_layers;
  // Cell store index of each cell of a layer, in the order of m_layer_cell_ids
  std::map<std::size_t, std::vector<std::uint32_t>> layer_store_index;

  // Iterate through all cells to rebuild the layer-to-cell mapping
  for (size_t i = 0; i < m_cell_store.size(); ++i) {
    const Cell cell = m_cell_store.get_at_index(i);
    if (!cell.is_valid()) {
      throw std::runtime_error("Invalid cell found in cell store");
    }
    std::size_t const layer_id = cell.layer();
    reject_r_phi_z_cell(cell.isRPhiZ(), layer_id);
    unique_layers.insert(layer_id);
    m_layer_cell_ids[layer_id].push_back(cell.id());
    layer_store_index[layer_id].push_back(static_cast<std::uint32_t>(i));
    m_n_total_cells++;

    // Rebuild layer flags if not already done
    if (!m_layer_flags.contains(layer_id)) {
      m_layer_flags[layer_id] = LayerFlags {};
      m_layer_flags[layer_id].is_barrel = cell.isBarrel();
      layer_lookup(layer_id).coordinate_system = coordinate_system_of(
          cell.isXYZ(), cell.isEtaPhiR(), cell.isEtaPhiZ());

      // Initialize eta extensions
      m_layer_flags[layer_id].eta_extensions[kEtaPositive].eta_min =
          std::numeric_limits<double>::max();
      m_layer_flags[layer_id].eta_extensions[kEtaPositive].eta_max =
          -std::numeric_limits<double>::max();
      m_layer_flags[layer_id].eta_extensions[kEtaNegative].eta_min =
          std::numeric_limits<double>::max();
      m_layer_flags[layer_id].eta_extensions[kEtaNegative].eta_max =
          -std::numeric_limits<double>::max();
    }

    // Update eta extremes for each cell
    update_eta_extremes(layer_id, cell);
  }

  m_layer_ids = sorted_layer_ids(unique_layers);

  // Build the RTree of each layer from its cells
  for (const auto& [layer_id, store_index] : layer_store_index) {
    std::vector<Cell> cells;
    cells.reserve(store_index.size());
    for (const auto index : store_index) {
      cells.push_back(m_cell_store.get_at_index(index));
    }

    const auto coord_sys = coordinate_system(layer_id);
    LayerLookup& lookup = layer_lookup(layer_id);
    lookup.rtree = std::make_unique<RTree>(coord_sys, cells);
    lookup.store_index = store_index;
  }

  auto end_time = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> const elapsed = end_time - start_time;

  // Print timing information
  // CaloGeo has no output level
  FCS_MSG_LVL_NOCHK(::FastCaloSim::MSG::INFO,
                    "Done loading calo geometry. Took "
                        << elapsed.count() << " s (" << m_n_total_cells
                        << " cells in " << n_layers() << " layers)");
}
