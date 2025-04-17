// Copyright (c) 2024 CERN for the benefit of the FastCaloSim project

#include "FastCaloSim/Geometry/CaloGeo.h"

CaloGeo::CaloGeo(ROOT::RDataFrame& geo)
{
  build(geo);
}

auto CaloGeo::get_cell(unsigned int layer, const Position& pos) const -> Cell
{
  // Check if an alternative geometry handler is set for the layer
  if (m_alt_geo_handlers.count(layer) != 0U) {
    auto cell_id = m_alt_geo_handlers.at(layer)->get_cell_id(layer, pos);

    // Return an invalid cell if the alternative geometry handler returns -1
    if (cell_id == -1) {
      Cell invalid_cell;
      return invalid_cell;
    }
    return m_cell_id_map.at(cell_id);
  }

  // Else proceed with the default geometry handler
  return m_layer_tree_map.at(layer).query_point(pos);
}

auto CaloGeo::get_cell_id(unsigned int layer,
                          const Position& pos) const -> unsigned long long
{
  auto cell = get_cell(layer, pos);
  return cell.id();
}

auto CaloGeo::n_cells(unsigned int layer) const -> unsigned int
{
  return m_layer_tree_map.at(layer).size();
}

auto CaloGeo::n_cells() const -> unsigned int
{
  return m_n_total_cells;
}

auto CaloGeo::n_layers() const -> unsigned int
{
  return m_n_layers;
}

auto CaloGeo::get_cell(unsigned long long id) const -> Cell
{
  return m_cell_id_map.at(id);
}

auto CaloGeo::get_cell_at_idx(unsigned int layer, size_t idx) -> Cell
{
  return m_layer_tree_map.at(layer).at(idx);
}

auto CaloGeo::is_barrel(unsigned int layer) const -> bool
{
  return m_layer_flags.at(layer).is_barrel;
}

auto CaloGeo::is_xyz(unsigned int layer) const -> bool
{
  return m_layer_flags.at(layer).is_xyz;
}

auto CaloGeo::is_eta_phi_r(unsigned int layer) const -> bool
{
  return m_layer_flags.at(layer).is_eta_phi_r;
}

auto CaloGeo::is_eta_phi_z(unsigned int layer) const -> bool
{
  return m_layer_flags.at(layer).is_eta_phi_z;
}

auto CaloGeo::zpos(unsigned int layer,
                   const Position& pos,
                   Cell::SubPos subpos) const -> double
{
  const auto& cell = get_cell(layer, pos);
  return cell.z(subpos);
}

auto CaloGeo::rpos(unsigned int layer,
                   const Position& pos,
                   Cell::SubPos subpos) const -> double
{
  const auto& cell = get_cell(layer, pos);
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
  for (int i = ilayer; i <= llayer; ++i) {
    set_alt_geo_handler(i, handler);
  }
}

void CaloGeo::record_cell(const Cell& cell)
{
  auto layer = cell.layer();
  m_layer_tree_map[layer].insert_cell(cell);
  m_cell_id_map.emplace(cell.id(), cell);
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

void CaloGeo::build(ROOT::RDataFrame& geo)
{
  // Layer information of cells
  auto layer = geo.Take<unsigned int>("layer");

  // Number of layers from unique layer values
  m_n_layers = std::set<unsigned int>(layer->begin(), layer->end()).size();

  // Initialize layer flags
  for (int i = 0; i < m_n_layers; ++i) {
    m_layer_flags.emplace(i, LayerFlags {});
  }
  // Count the total number of cells
  m_n_total_cells = *geo.Count();

  auto isBarrel = geo.Take<bool>("isBarrel");
  auto id = geo.Take<unsigned long long>("id");
  auto x = geo.Take<float>("x");
  auto y = geo.Take<float>("y");
  auto z = geo.Take<float>("z");
  auto phi = geo.Take<float>("phi");
  auto eta = geo.Take<float>("eta");
  auto r = geo.Take<float>("r");
  auto dx = geo.Take<float>("dx");
  auto dy = geo.Take<float>("dy");
  auto dz = geo.Take<float>("dz");
  auto dphi = geo.Take<float>("dphi");
  auto deta = geo.Take<float>("deta");
  auto dr = geo.Take<float>("dr");
  auto isXYZ = geo.Take<bool>("isXYZ");
  auto isEtaPhiR = geo.Take<bool>("isEtaPhiR");
  auto isEtaPhiZ = geo.Take<bool>("isEtaPhiZ");
  auto isRPhiZ = geo.Take<bool>("isRPhiZ");

  for (size_t i = 0; i < m_n_total_cells; ++i) {
    m_layer_flags.at(layer->at(i)).is_barrel = isBarrel->at(i);
    m_layer_flags.at(layer->at(i)).is_xyz = isXYZ->at(i);
    m_layer_flags.at(layer->at(i)).is_eta_phi_r = isEtaPhiR->at(i);
    m_layer_flags.at(layer->at(i)).is_eta_phi_z = isEtaPhiZ->at(i);
    m_layer_flags.at(layer->at(i)).is_r_phi_z = isRPhiZ->at(i);

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
    if (isRPhiZ->at(i) && deta->at(i) <= 0.0) {
      std::array<float, 4> etas;
      int corner_idx = 0;
      for (float r_sign : {-0.5, 0.5}) {
        for (float z_sign : {-0.5, 0.5}) {
          float r_val = r->at(i) + r_sign * dr->at(i);
          float z_val = z->at(i) + z_sign * dz->at(i);
          float theta = std::atan2(r_val, z_val);
          etas[corner_idx++] = -std::log(std::tan(theta / 2));
        }
      }

      auto minmax = std::minmax_element(etas.begin(), etas.end());
      deta->at(i) = *minmax.second - *minmax.first;
    }

    Position pos = Position {.m_x = x->at(i),
                             .m_y = y->at(i),
                             .m_z = z->at(i),
                             .m_eta = eta->at(i),
                             .m_phi = phi->at(i),
                             .m_r = r->at(i)};

    Cell cell = Cell {id->at(i),
                      pos,
                      layer->at(i),
                      isBarrel->at(i),
                      isXYZ->at(i),
                      isEtaPhiR->at(i),
                      isEtaPhiZ->at(i),
                      isRPhiZ->at(i),
                      dx->at(i),
                      dy->at(i),
                      dz->at(i),
                      deta->at(i),
                      dphi->at(i),
                      dr->at(i)};

    record_cell(cell);
    update_eta_extremes(layer->at(i), cell);
  }
}
