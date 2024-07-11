#include <map>
#include <vector>

#include <ROOT/RDataFrame.hxx>

#include "FastCaloSim/Geometry/Cell.h"
#include "FastCaloSim/Geometry/RTree3D.h"

class CaloCellRDTree
{
public:
  CaloCellRDTree(ROOT::RDataFrame& geo) { build(geo); }

  auto get_num_cells(int layer) const { return m_tree_map.at(layer).size(); }

  auto get_cell(int layer, double x, double y, double z) -> Cell
  {
    return m_tree_map.at(layer).query_point(x, y, z);
  }

  auto get_cell_at_index(int layer, size_t idx) const
  {
    return m_tree_map.at(layer).at(idx);
  }

private:
  void build(ROOT::RDataFrame& geo)
  {
    size_t numEntries = *geo.Count();

    std::map<int, std::vector<Cell>> layer_cells_map;

    auto layer_data = geo.Take<long long>("layer");
    auto id_data = geo.Take<long long>("id");
    auto x_data = geo.Take<double>("x");
    auto y_data = geo.Take<double>("y");
    auto z_data = geo.Take<double>("z");
    auto phi_data = geo.Take<double>("phi");
    auto eta_data = geo.Take<double>("eta");
    auto r_data = geo.Take<double>("r");
    auto dx_data = geo.Take<double>("dx");
    auto dy_data = geo.Take<double>("dy");
    auto dz_data = geo.Take<double>("dz");
    auto dphi_data = geo.Take<double>("dphi");
    auto deta_data = geo.Take<double>("deta");
    auto dr_data = geo.Take<double>("dr");
    auto isBarrel_data = geo.Take<long long>("isBarrel");
    auto isCartesian_data = geo.Take<long long>("isCartesian");
    auto isCylindrical_data = geo.Take<long long>("isCylindrical");
    auto isECCylindrical_data = geo.Take<long long>("isECCylindrical");

    for (size_t i = 0; i < numEntries; ++i) {
      Cell cell = {id_data->at(i),
                   x_data->at(i),
                   y_data->at(i),
                   z_data->at(i),
                   phi_data->at(i),
                   eta_data->at(i),
                   r_data->at(i),
                   dx_data->at(i),
                   dy_data->at(i),
                   dz_data->at(i),
                   dphi_data->at(i),
                   deta_data->at(i),
                   dr_data->at(i),
                   static_cast<bool>(isBarrel_data->at(i)),
                   static_cast<bool>(isCartesian_data->at(i)),
                   static_cast<bool>(isCylindrical_data->at(i)),
                   static_cast<bool>(isECCylindrical_data->at(i)),
                   static_cast<int>(layer_data->at(i))};

      // Insert cell into the layer map
      layer_cells_map[static_cast<int>(layer_data->at(i))].push_back(cell);
    }

    for (const auto& [layer, cells] : layer_cells_map) {
      std::cout << "Building RTree for layer... " << layer << '\n';
      for (const auto& cell : cells) {
        m_tree_map[layer].insert_cell(cell);
      }
    }
  }

  std::map<int, RTree3D> m_tree_map;
};
