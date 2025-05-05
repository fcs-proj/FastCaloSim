// Copyright (c) 2025 CERN for the benefit of the FastCaloSim project

#include <array>
#include <chrono>
#include <cmath>
#include <iostream>
#include <limits>
#include <optional>
#include <unordered_set>
#include <vector>

#include "FastCaloSim/Geometry/CaloGeo.h"

#include "FastCaloSim/Geometry/CellStoreBuilder.h"

auto CaloGeo::get_cell(unsigned int layer,
                       const Position& pos) const -> const Cell&
{
  // Check if an alternative geometry handler is set for the layer
  auto alt_it = m_alt_geo_handlers.find(layer);
  if (alt_it != m_alt_geo_handlers.end()) {
    auto cell_id = alt_it->second->get_cell_id(layer, pos);
    // Return an invalid cell if the alternative geometry handler returns -1
    if (cell_id == std::numeric_limits<unsigned long long>::max()) {
      std::runtime_error("Invalid cell ID from alternative geometry handler");
      static Cell invalid_cell;
      return invalid_cell;
    }
    return get_cell(cell_id);
  }

  // Else proceed with the default geometry handler
  auto query_it = m_layer_rtree_queries.find(layer);
  if (query_it == m_layer_rtree_queries.end()) {
    throw std::runtime_error("No RTree loaded for layer "
                             + std::to_string(layer));
  }

  // Query the RTree to get the cell ID directly
  auto cell_id = query_it->second->query_point(pos);

  // Look up the cell in the repository
  return get_cell(cell_id);
}

auto CaloGeo::get_cell_id(unsigned int layer,
                          const Position& pos) const -> unsigned long long
{
  const Cell& cell = get_cell(layer, pos);
  return cell.id();
}

auto CaloGeo::n_cells() const -> unsigned int
{
  return m_n_total_cells;
}

auto CaloGeo::n_cells(unsigned int layer) const -> unsigned int
{
  auto it = m_layer_cell_ids.find(layer);
  if (it != m_layer_cell_ids.end()) {
    return static_cast<unsigned int>(it->second.size());
  }
  return 0;
}

auto CaloGeo::n_layers() const -> unsigned int
{
  return m_n_layers;
}

auto CaloGeo::get_cell(unsigned long long id) const -> const Cell&
{
  return m_cell_store.get(id);
}

auto CaloGeo::get_cell_at_idx(unsigned int layer,
                              size_t idx) const -> const Cell&
{
  auto layer_it = m_layer_cell_ids.find(layer);
  if (layer_it == m_layer_cell_ids.end() || idx >= layer_it->second.size()) {
    throw std::runtime_error("Invalid layer ID or index out of bounds");
  }

  unsigned long long cell_id = layer_it->second[idx];
  return m_cell_store.get(cell_id);
}

auto CaloGeo::is_barrel(unsigned int layer) const -> bool
{
  return m_layer_flags.at(layer).is_barrel;
}

auto CaloGeo::is_xyz(unsigned int layer) const -> bool
{
  return m_layer_flags.at(layer).coordinate_system
      == RTreeHelpers::CoordinateSystem::XYZ;
}

auto CaloGeo::is_eta_phi_r(unsigned int layer) const -> bool
{
  return m_layer_flags.at(layer).coordinate_system
      == RTreeHelpers::CoordinateSystem::EtaPhiR;
}

auto CaloGeo::is_eta_phi_z(unsigned int layer) const -> bool
{
  return m_layer_flags.at(layer).coordinate_system
      == RTreeHelpers::CoordinateSystem::EtaPhiZ;
}

auto CaloGeo::zpos(unsigned int layer,
                   const Position& pos,
                   Cell::SubPos subpos) const -> double
{
  const Cell& cell = get_cell(layer, pos);
  return cell.z(subpos);
}

auto CaloGeo::rpos(unsigned int layer,
                   const Position& pos,
                   Cell::SubPos subpos) const -> double
{
  const Cell& cell = get_cell(layer, pos);
  return cell.r(subpos);
}

auto CaloGeo::min_eta(unsigned int layer, DetectorSide side) const -> double
{
  return m_layer_flags.at(layer).eta_extensions.at(side).eta_min;
}

auto CaloGeo::max_eta(unsigned int layer, DetectorSide side) const -> double
{
  return m_layer_flags.at(layer).eta_extensions.at(side).eta_max;
}

void CaloGeo::set_alt_geo_handler(unsigned int layer,
                                  std::shared_ptr<CaloGeo> handler)
{
  m_alt_geo_handlers[layer] = std::move(handler);
}

void CaloGeo::set_alt_geo_handler(unsigned int ilayer,
                                  unsigned int llayer,
                                  std::shared_ptr<CaloGeo> handler)
{
  for (unsigned int i = ilayer; i <= llayer; ++i) {
    m_alt_geo_handlers[i] = handler;
  }
}

void CaloGeo::update_eta_extremes(unsigned int layer, const Cell& cell)
{
  DetectorSide side = cell.eta() > 0 ? kEtaPositive : kEtaNegative;
  double half_width = cell.deta() / 2;
  double min_eta = cell.eta() - half_width;
  double max_eta = cell.eta() + half_width;
  auto& eta_extremes = m_layer_flags.at(layer).eta_extensions.at(side);
  if (max_eta > eta_extremes.eta_max) {
    eta_extremes.eta_max = max_eta;
  }
  if (min_eta < eta_extremes.eta_min) {
    eta_extremes.eta_min = min_eta;
  }
}

// Method to build geometry
void CaloGeo::build(ROOT::RDataFrame& geo, const std::string& rtree_base_path)
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
  std::unordered_set<unsigned int> unique_layers;
  for (const auto& layer_id : *layer) {
    unique_layers.insert(layer_id);
  }
  m_n_layers = unique_layers.size();

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
  std::unordered_set<unsigned int> processed_layer_flags;
  // This will be used to build our memory-mapped cell store
  CellStoreBuilder cell_store_builder;

  // Process cells in bulk
  Position pos;
  for (size_t i = 0; i < m_n_total_cells; ++i) {
    unsigned int layer_id = layer->at(i);

    // Only update layer flags once per layer
    if (processed_layer_flags.find(layer_id) == processed_layer_flags.end()) {
      auto& flags = m_layer_flags[layer_id];
      flags.is_barrel = isBarrel->at(i);
      // Set coordinate system based on the boolean flags
      if (isXYZ->at(i)) {
        flags.coordinate_system = RTreeHelpers::CoordinateSystem::XYZ;
      } else if (isEtaPhiR->at(i)) {
        flags.coordinate_system = RTreeHelpers::CoordinateSystem::EtaPhiR;
      } else if (isEtaPhiZ->at(i)) {
        flags.coordinate_system = RTreeHelpers::CoordinateSystem::EtaPhiZ;
      } else if (isRPhiZ->at(i)) {
        flags.coordinate_system = RTreeHelpers::CoordinateSystem::RPhiZ;
      } else {
        flags.coordinate_system = RTreeHelpers::CoordinateSystem::Undefined;
      }
      processed_layer_flags.insert(layer_id);
    }

    // For RPhiZ cells, we estimate the cell's Δη if it is missing (≤ 0).
    // RPhiZ cells are defined in cylindrical coordinates (r, φ, z)
    // and do not have η as a primary coordinate. If deta is not provided in the
    // input, we compute an approximate Δη by evaluating η at the four (r, z)
    // corners of the cell and taking the span (η_max - η_min). This gives a
    // good estimate of the angular extent, especially for projective or
    // wedge-like cell shapes where both r and z vary.
    //
    // We do not apply this approximation to EtaPhiZ or EtaPhiR cells,
    // since those coordinate systems use η as a fundamental axis
    // and should have accurate Δη values already present in the input.
    float deta_val = deta->at(i);
    if (isRPhiZ->at(i) && deta_val <= 0.0f) {
      std::array<double, 4> etas;
      for (int corner = 0; corner < 4; ++corner) {
        double r_sign = (corner & 1) ? 0.5 : -0.5;
        double z_sign = (corner & 2) ? 0.5 : -0.5;

        double r_val = r->at(i) + r_sign * dr->at(i);
        double z_val = z->at(i) + z_sign * dz->at(i);
        double theta = std::atan2(r_val, z_val);
        etas[corner] = -std::log(std::tan(theta / 2));
      }

      auto minmax = std::minmax_element(etas.begin(), etas.end());
      deta_val = static_cast<float>(*minmax.second - *minmax.first);
    }

    pos.m_x = x->at(i);
    pos.m_y = y->at(i);
    pos.m_z = z->at(i);
    pos.m_eta = eta->at(i);
    pos.m_phi = phi->at(i);
    pos.m_r = r->at(i);

    Cell cell(id->at(i),
              pos,
              layer_id,
              isBarrel->at(i),
              isXYZ->at(i),
              isEtaPhiR->at(i),
              isEtaPhiZ->at(i),
              isRPhiZ->at(i),
              dx->at(i),
              dy->at(i),
              dz->at(i),
              deta_val,
              dphi->at(i),
              dr->at(i));

    update_eta_extremes(layer_id, cell);
    m_layer_cell_ids[layer_id].push_back(cell.id());

    cell_store_builder.add_cell(cell);
  }

  // Write the cell store to disk
  // Memory will be freed after writing
  cell_store_builder.write(rtree_base_path + "/cellstore");

  // Temporary cell store
  CellStore temp_cell_store;
  temp_cell_store.load(rtree_base_path + "/cellstore");

  // Build RTrees for each layer
  for (const auto& [layer_id, cell_ids] : m_layer_cell_ids) {
    // Get coordinate system directly from layer flags
    auto coord_sys = m_layer_flags.at(layer_id).coordinate_system;

    std::string rtree_name = "rtree_layer_" + std::to_string(layer_id);
    std::string rtree_path = rtree_base_path + "/" + rtree_name;

    // Build the RTree
    RTreeBuilder rtree_builder(coord_sys);

    for (size_t i = 0; i < cell_ids.size(); ++i) {
      const Cell& cell = temp_cell_store.get(cell_ids[i]);
      rtree_builder.add_cell(&cell);
    }
    rtree_builder.build(rtree_path);

    std::cout << "Built RTree for layer " << layer_id << " with "
              << cell_ids.size() << " cells" << std::endl;
  }

  auto end_time = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = end_time - start_time;

  // Print timing information
  std::cout << "INFO: Done building calo geometry. Took " << elapsed.count()
            << " s (" << m_n_total_cells << " cells in " << m_n_layers
            << " layers)" << std::endl;
}

// Method to load geometry
void CaloGeo::load(
    const std::string& rtree_base_path,
    size_t rtree_cache_size,  // per-layer r-tree cache size in bytes
    size_t cell_store_cache_size)  // cell store cache size in bytes
{
  // Start timing
  auto start_time = std::chrono::high_resolution_clock::now();

  // Load the cell store
  m_cell_store.load(rtree_base_path + "/cellstore", cell_store_cache_size);

  // Rebuild the layer cell ID map
  m_layer_cell_ids.clear();
  m_n_total_cells = 0;
  std::unordered_set<unsigned int> unique_layers;

  // Iterate through all cells to rebuild the layer-to-cell mapping
  for (size_t i = 0; i < m_cell_store.size(); ++i) {
    const Cell& cell = m_cell_store.get_at_index(i);
    if (!cell.is_valid()) {
      throw std::runtime_error("Invalid cell found in cell store");
    }
    unsigned int layer_id = cell.layer();
    unique_layers.insert(layer_id);
    m_layer_cell_ids[layer_id].push_back(cell.id());
    m_n_total_cells++;

    // Rebuild layer flags if not already done
    if (m_layer_flags.find(layer_id) == m_layer_flags.end()) {
      m_layer_flags[layer_id] = LayerFlags {};
      m_layer_flags[layer_id].is_barrel = cell.isBarrel();

      // Set coordinate system based on cell flags
      if (cell.isXYZ()) {
        m_layer_flags[layer_id].coordinate_system =
            RTreeHelpers::CoordinateSystem::XYZ;
      } else if (cell.isEtaPhiR()) {
        m_layer_flags[layer_id].coordinate_system =
            RTreeHelpers::CoordinateSystem::EtaPhiR;
      } else if (cell.isEtaPhiZ()) {
        m_layer_flags[layer_id].coordinate_system =
            RTreeHelpers::CoordinateSystem::EtaPhiZ;
      } else if (cell.isRPhiZ()) {
        m_layer_flags[layer_id].coordinate_system =
            RTreeHelpers::CoordinateSystem::RPhiZ;
      } else {
        m_layer_flags[layer_id].coordinate_system =
            RTreeHelpers::CoordinateSystem::Undefined;
      }

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

  m_n_layers = unique_layers.size();

  // Load the RTrees for each layer
  for (const auto& [layer_id, cell_ids] : m_layer_cell_ids) {
    auto coord_sys = m_layer_flags.at(layer_id).coordinate_system;

    std::string rtree_name = "rtree_layer_" + std::to_string(layer_id);
    std::string rtree_path = rtree_base_path + "/" + rtree_name;

    // Load the RTree for querying
    m_layer_rtree_queries[layer_id] = std::make_unique<RTreeQuery>(coord_sys);
    // Load the RTree from disk with specified cache size
    m_layer_rtree_queries[layer_id]->load(rtree_path, rtree_cache_size);
  }

  auto end_time = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = end_time - start_time;

  // Print timing information
  std::cout << "INFO: Done loading calo geometry. Took " << elapsed.count()
            << " s (" << m_n_total_cells << " cells in " << m_n_layers
            << " layers)" << std::endl;
}
