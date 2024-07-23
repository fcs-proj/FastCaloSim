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
    m_tree_map[cell->layer()].insert_cell(cell);
  }

  void build(ROOT::RDataFrame& geo)
  {
    auto numEntries = *geo.Count();

    auto layer = geo.Take<long long>("layer");
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

    // Loop over cells and build RTree
    for (size_t i = 0; i < numEntries; ++i) {
      if (isCartesian->at(i)) {
        XYZCell cell = {id->at(i),
                        x->at(i),
                        y->at(i),
                        z->at(i),
                        phi->at(i),
                        eta->at(i),
                        r->at(i),
                        layer->at(i),
                        dx->at(i),
                        dy->at(i),
                        dz->at(i)};

        record_cell(cell);

      } else if (isCylindrical->at(i)) {
        EtaPhiRCell cell = {id->at(i),
                            x->at(i),
                            y->at(i),
                            z->at(i),
                            phi->at(i),
                            eta->at(i),
                            r->at(i),
                            layer->at(i),
                            deta->at(i),
                            dphi->at(i),
                            dr->at(i)};

        record_cell(cell);
      } else if (isECCylindrical->at(i)) {
        EtaPhiZCell cell = {id->at(i),
                            x->at(i),
                            y->at(i),
                            z->at(i),
                            phi->at(i),
                            eta->at(i),
                            r->at(i),
                            layer->at(i),
                            deta->at(i),
                            dphi->at(i),
                            dz->at(i)};
        record_cell(cell);
      }
    }
  }

  std::map<int, RTree> m_tree_map;
};
