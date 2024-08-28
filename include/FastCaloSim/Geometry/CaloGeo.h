#pragma once

#include <map>
#include <vector>

#include <ROOT/RDataFrame.hxx>

#include "FastCaloSim/Geometry/Cell.h"
#include "FastCaloSim/Geometry/RTree.h"

class CaloGeo
{
public:
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

  // @brief Returns the best matching cell for a hit in a specific layer
  // Note that the hit must have the same coordinate system as the cells in the
  // tree
  template<typename T>
  auto get_cell(int layer, const T& hit) const -> Cell
  {
    return m_layer_tree_map.at(layer).query_point(hit);
  }

  // @brief Returns a boolean indicating if a layer is a barrel or endcap
  auto is_barrel(int layer) const { return m_is_layer_barrel_map.at(layer); }

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

  // @brief returns the r of the cell with maximum r
  // For endcap layers, subpos will be ignored and the center position returned
  auto rmax(int layer, Cell::SubPos subpos) const -> double
  {
    auto cell = get_cell(m_edge_cells_ids.at(layer).rmax);
    return cell.r(subpos);
  }
  // @brief returns the r of the cell with minimum r
  // For endcap layers, subpos will be ignored and the center position returned
  auto rmin(int layer, Cell::SubPos subpos) const -> double
  {
    auto cell = get_cell(m_edge_cells_ids.at(layer).rmin);
    return cell.r(subpos);
  }

  // @brief returns the z of the cell with maximum z
  // For barrel layers, subpos will be ignored and the center position returned
  auto zmax(int layer, Cell::SubPos subpos) const -> double
  {
    auto cell = get_cell(m_edge_cells_ids.at(layer).zmax);
    return cell.z(subpos);
  }
  // @brief returns the z of the cell with minimum z
  // For barrel layers, subpos will be ignored and the center position returned
  auto zmin(int layer, Cell::SubPos subpos) const -> double
  {
    auto cell = get_cell(m_edge_cells_ids.at(layer).zmin);
    return cell.z(subpos);
  }
  // @brief returns the eta of the edge of the cell cell with minimum eta in a
  // layer For cells in xzy coordinates, this will return the center position
  auto etamin(int layer) const -> double
  {
    auto cell = get_cell(m_edge_cells_ids.at(layer).etamin);

    if (cell.eta() > 0)
      std::runtime_error("Minimum cell eta should not be positive");

    return cell.isXYZ() ? cell.eta() : cell.eta() + cell.deta() * 0.5;
  }

  // @brief returns the eta of the cell with maximum eta in a layer
  // For cells in xzy coordinates, this will return the center position
  auto etamax(int layer) const -> double
  {
    auto cell = get_cell(m_edge_cells_ids.at(layer).etamax);

    if (cell.eta() < 0)
      std::runtime_error("Maximum cell eta should not be negative");

    return cell.isXYZ() ? cell.eta() : cell.eta() - cell.deta() * 0.5;
  }

private:
  struct EdgeCellIds
  {
    long long rmin, rmax;
    long long zmin, zmax;
    long long etamin, etamax;
  };

  int m_n_layers {};
  int m_n_total_cells {};
  std::map<int, RTree> m_layer_tree_map;
  std::map<long long, Cell> m_cell_id_map;
  std::unordered_map<int, bool> m_is_layer_barrel_map;
  std::unordered_map<int, EdgeCellIds> m_edge_cells_ids;

  // @brief Record the cell in the geometry
  void record_cell(const Cell& cell)
  {
    m_layer_tree_map[cell.layer()].insert_cell(cell);
    m_cell_id_map.emplace(cell.id(), cell);
  }

  // @brief Update the edge cells for a layer
  void update_edge_cells(int layer, const Cell& cell)
  {
    // Initialize for first cell in a layer
    if (m_edge_cells_ids.find(layer) == m_edge_cells_ids.end()) {
      m_edge_cells_ids[layer] = EdgeCellIds {
          cell.id(), cell.id(), cell.id(), cell.id(), cell.id(), cell.id()};
    }
    auto& ext = m_edge_cells_ids.at(layer);

    const auto& rmax_cell = get_cell(ext.rmax);
    const auto& rmin_cell = get_cell(ext.rmin);
    const auto& zmax_cell = get_cell(ext.zmax);
    const auto& zmin_cell = get_cell(ext.zmin);
    const auto& etamax_cell = get_cell(ext.etamax);
    const auto& etamin_cell = get_cell(ext.etamin);

    ext.rmax = cell.r() > rmax_cell.r() ? cell.id() : ext.rmax;
    ext.rmin = cell.r() < rmin_cell.r() ? cell.id() : ext.rmin;
    ext.zmax = cell.z() > zmax_cell.z() ? cell.id() : ext.zmax;
    ext.zmin = cell.z() < zmin_cell.z() ? cell.id() : ext.zmin;
    ext.etamax = cell.eta() > etamax_cell.eta() ? cell.id() : ext.etamax;
    ext.etamin = cell.eta() < etamin_cell.eta() ? cell.id() : ext.etamin;
  }

  // @brief Build the geometry from the RDataFrame
  void build(ROOT::RDataFrame& geo)
  {
    // Record number of layers
    m_n_layers = geo.Max<long long>("layer").GetValue();

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

      m_is_layer_barrel_map.emplace(layer->at(i), is_barrel);

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
      update_edge_cells(layer->at(i), cell);
    }
  }
};
