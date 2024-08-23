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
  auto n_cells() const { return m_n_cells; }

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
  auto is_barrel(int layer) const { return m_layer_flags.at(layer).is_barrel; }

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
  struct LayerFlags
  {
    bool is_barrel;
    bool is_xyz, is_EtaPhiR, is_EtaPhiZ;
  };

  struct EdgeCellIds
  {
    long long rmin, rmax;
    long long zmin, zmax;
    long long etamin, etamax;
  };

  int m_n_layers {};
  int m_n_cells {};
  std::map<int, RTree> m_layer_tree_map;
  std::map<long long, Cell> m_cell_id_map;
  std::unordered_map<int, LayerFlags> m_layer_flags;
  std::unordered_map<int, EdgeCellIds> m_edge_cells_ids;

  // @Brief Returns the edge cells of a layer with maximum and minimum r and z
  static auto record_edge_cell_ids(ROOT::RDataFrame& df, int layer)
      -> EdgeCellIds
  {
    // Filter the DataFrame to include only the specified layer
    auto layer_df =
        df.Filter([layer](long long l) { return l == layer; }, {"layer"});

    // Lambda to get the id of the cell with a specific extreme value in a
    // column
    auto get_id_for_extreme = [&](const std::string& col, double extreme_value)
    {
      return layer_df
          .Filter([extreme_value](double val) { return val == extreme_value; },
                  {col})
          .Take<long long>("id")
          ->front();
    };

    // Find the min and max values for r and z
    auto [rmin, rmax] = std::make_tuple(*layer_df.Min("r"), *layer_df.Max("r"));
    auto [zmin, zmax] = std::make_tuple(*layer_df.Min("z"), *layer_df.Max("z"));
    auto [etamin, etamax] =
        std::make_tuple(*layer_df.Min("eta"), *layer_df.Max("eta"));

    // Get the corresponding ids for the min/max values
    return {get_id_for_extreme("r", rmin),
            get_id_for_extreme("r", rmax),
            get_id_for_extreme("z", zmin),
            get_id_for_extreme("z", zmax),
            get_id_for_extreme("eta", etamin),
            get_id_for_extreme("eta", etamax)};
  }

  void record_cell(const Cell& cell)
  {
    m_layer_tree_map[cell.layer()].insert_cell(cell);
    m_cell_id_map.emplace(cell.id(), cell);
  }

  void build(ROOT::RDataFrame& geo)
  {
    std::cout << "Building RTree geometry..." << '\n';

    // Record number of layers
    m_n_layers = geo.Max<long long>("layer").GetValue();

    // Record number of cells
    m_n_cells = *geo.Count();

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

    for (size_t i = 0; i < m_n_cells; ++i) {
      m_layer_flags.emplace(
          layer->at(i),
          LayerFlags {static_cast<bool>(isBarrel->at(i)),
                      static_cast<bool>(isCartesian->at(i)),
                      static_cast<bool>(isCylindrical->at(i)),
                      static_cast<bool>(isECCylindrical->at(i))});

      // Position of the center of the cell
      Position pos = Position {.m_x = x->at(i),
                               .m_y = y->at(i),
                               .m_z = z->at(i),
                               .m_eta = eta->at(i),
                               .m_phi = phi->at(i),
                               .m_r = r->at(i)};

      record_cell(Cell {id->at(i),
                        pos,
                        layer->at(i),
                        m_layer_flags.at(layer->at(i)).is_barrel,
                        m_layer_flags.at(layer->at(i)).is_xyz,
                        m_layer_flags.at(layer->at(i)).is_EtaPhiR,
                        m_layer_flags.at(layer->at(i)).is_EtaPhiZ,
                        dx->at(i),
                        dy->at(i),
                        dz->at(i),
                        deta->at(i),
                        dphi->at(i),
                        dr->at(i)});
    }

    // Record the cells at the edge of each layer (with maximum r and z center
    // positions)
    for (int layer = 0; layer <= m_n_layers; ++layer) {
      m_edge_cells_ids[layer] = record_edge_cell_ids(geo, layer);
    }

    std::cout << "Building geometry completed!" << '\n';
  }
};
