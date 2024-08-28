#pragma once

#include <map>
#include <vector>

#include <ROOT/RDataFrame.hxx>

#include "FastCaloSim/Geometry/Cell.h"
#include "FastCaloSim/Geometry/RTree.h"

class CaloCellLookup
{
public:
  CaloCellLookup(ROOT::RDataFrame& geo) { build(geo); }

  auto get_num_cells(int layer) { return m_tree_map.at(layer).size(); }

  template<typename T>
  auto get_cell(int layer, T& hit) -> Cell
  {
    return m_tree_map.at(layer).query_point(hit);
  }

  auto get_cell_at_index(int layer, size_t idx) const
  {
    return m_tree_map.at(layer).at(idx);
  }

private:
  void record_cell(const Cell& cell)
  {
    m_tree_map[cell.layer()].insert_cell(cell);
  }

  void build(ROOT::RDataFrame& geo)
  {
    auto numEntries = *geo.Count();

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

    for (size_t i = 0; i < numEntries; ++i) {
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
                        static_cast<bool>(isBarrel->at(i)),
                        static_cast<bool>(isCartesian->at(i)),
                        static_cast<bool>(isCylindrical->at(i)),
                        static_cast<bool>(isECCylindrical->at(i)),
                        dx->at(i),
                        dy->at(i),
                        dz->at(i),
                        deta->at(i),
                        dphi->at(i),
                        dr->at(i)});
    }
  }

  std::map<int, RTree> m_tree_map;
};
