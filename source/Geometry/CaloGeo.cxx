// Copyright (c) 2024 CERN for the benefit of the FastCaloSim project

#include "FastCaloSim/Geometry/CaloGeo.h"

CaloGeo::CaloGeo(ROOT::RDataFrame& geo)
{
  build(geo);
}

auto CaloGeo::get_cell(int layer, const Position& pos) const -> Cell
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

auto CaloGeo::get_cell_id(int layer, const Position& pos) const -> long long
{
  auto cell = get_cell(layer, pos);
  return cell.id();
}

auto CaloGeo::n_cells(int layer) const -> unsigned int
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

auto CaloGeo::get_cell(long long id) const -> Cell
{
  return m_cell_id_map.at(id);
}

auto CaloGeo::get_cell_at_idx(int layer, size_t idx) -> Cell
{
  return m_layer_tree_map.at(layer).at(idx);
}

auto CaloGeo::is_barrel(int layer) const -> bool
{
  return m_layer_flags.at(layer).is_barrel;
}

auto CaloGeo::is_xyz(int layer) const -> bool
{
  return m_layer_flags.at(layer).is_xyz;
}

auto CaloGeo::is_eta_phi_r(int layer) const -> bool
{
  return m_layer_flags.at(layer).is_eta_phi_r;
}

auto CaloGeo::is_eta_phi_z(int layer) const -> bool
{
  return m_layer_flags.at(layer).is_eta_phi_z;
}

auto CaloGeo::zpos(int layer,
                   const Position& pos,
                   Cell::SubPos subpos) const -> double
{
  const auto& cell = get_cell(layer, pos);
  return cell.z(subpos);
}

auto CaloGeo::rpos(int layer,
                   const Position& pos,
                   Cell::SubPos subpos) const -> double
{
  const auto& cell = get_cell(layer, pos);
  return cell.r(subpos);
}

auto CaloGeo::min_eta(int layer, DetectorSide side) const -> double
{
  return m_layer_flags.at(layer).eta_extensions.at(side).eta_min;
}

auto CaloGeo::max_eta(int layer, DetectorSide side) const -> double
{
  return m_layer_flags.at(layer).eta_extensions.at(side).eta_max;
}

void CaloGeo::set_alt_geo_handler(int layer, std::shared_ptr<CaloGeo> handler)
{
  m_alt_geo_handlers[layer] = std::move(handler);
}

void CaloGeo::set_alt_geo_handler(int ilayer,
                                  int llayer,
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

void CaloGeo::update_eta_extremes(int layer, const Cell& cell)
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
  auto layer = geo.Take<long long>("layer");

  // Number of layers from unique layer values
  m_n_layers = std::set<long long>(layer->begin(), layer->end()).size();

  // Initialize layer flags
  for (int i = 0; i < m_n_layers; ++i) {
    m_layer_flags.emplace(i, LayerFlags {});
  }
  // Count the total number of cells
  m_n_total_cells = *geo.Count();

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

  for (size_t i = 0; i < m_n_total_cells; ++i) {
    bool is_barrel = static_cast<bool>(isBarrel->at(i));
    bool is_XYZ = static_cast<bool>(isCartesian->at(i));
    bool is_EtaPhiR = static_cast<bool>(isCylindrical->at(i));
    bool is_EtaPhiZ = static_cast<bool>(isECCylindrical->at(i));

    m_layer_flags.at(layer->at(i)).is_barrel = is_barrel;
    m_layer_flags.at(layer->at(i)).is_xyz = is_XYZ;
    m_layer_flags.at(layer->at(i)).is_eta_phi_r = is_EtaPhiR;
    m_layer_flags.at(layer->at(i)).is_eta_phi_z = is_EtaPhiZ;

    Position pos = Position {.m_x = x->at(i),
                             .m_y = y->at(i),
                             .m_z = z->at(i),
                             .m_eta = eta->at(i),
                             .m_phi = phi->at(i),
                             .m_r = r->at(i)};

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

    record_cell(cell);
    update_eta_extremes(layer->at(i), cell);
  }
}
