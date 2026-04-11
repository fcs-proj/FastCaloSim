// Copyright (c) 2025 CERN for the benefit of the FastCaloSim project
#include <memory>
#include <stdexcept>

#include "FastCaloSim/Geometry/RTreeBuilder.h"

using namespace FastCaloSim::Geometry;

RTreeBuilder::RTreeBuilder(RTreeHelpers::CoordinateSystem coordSys)
    : m_coordinate_system(coordSys)
{
}

void RTreeBuilder::add_cell(const Cell* cell)
{
  // Use the cell ID instead of the pointer
  uint64_t cell_id = cell->id();

  if (m_coordinate_system == RTreeHelpers::CoordinateSystem::XYZ) {
    auto box = RTreeHelpers::build_xy_box(
        cell->x(), cell->y(), cell->dx(), cell->dy());
    m_boxes.emplace_back(box, cell_id);
  } else {
    auto boxes = RTreeHelpers::build_eta_phi_boxes(
        cell->eta(), cell->phi(), cell->deta(), cell->dphi());
    for (const auto& box : boxes) {
      m_boxes.emplace_back(box, cell_id);
    }
  }
}

void RTreeBuilder::build(const std::string& output_path)
{
  if (m_boxes.empty()) {
    throw std::logic_error(
        "RTreeBuilder::build() called with no cells. "
        "Insert cells first via add_cell().");
  }

  // Create data stream for bulk loading
  class DataStream : public SpatialIndex::IDataStream
  {
  public:
    DataStream(
        const std::vector<std::pair<std::array<double, 4>, uint64_t>>& boxes)
        : m_boxes(boxes)
        , m_index(0)
    {
    }

    bool hasNext() override { return m_index < m_boxes.size(); }

    SpatialIndex::IData* getNext() override
    {
      // Get the next box-cell id pair from our collection
      const auto& [bounds, id] = m_boxes[m_index++];

      // Create a spatial region (bounding box) from the bounds
      double low[2] = {bounds[0], bounds[1]};  // Min coordinates
      double high[2] = {bounds[2], bounds[3]};  // Max coordinates
      SpatialIndex::Region r(low, high, 2);  // 2D region

      // Create a new Data object that associates the region with our cell
      // pointer The RTree will take ownership of this object and free it when
      // done
      return new SpatialIndex::RTree::Data(
          0,  // Size of our data (not used as we just return id
          nullptr,  // Pointer to data (not used as we just return id
          r,  // Region (bounding box)
          id  // ID for this entry
      );
    }

    void rewind() override
    {
      throw Tools::NotSupportedException("Rewind not supported");
    }

    uint32_t size() override { return static_cast<uint32_t>(m_boxes.size()); }

  private:
    const std::vector<std::pair<std::array<double, 4>, uint64_t>>& m_boxes;
    size_t m_index;
  };

  DataStream stream(m_boxes);

  // Create disk storage manager
  std::string filename = output_path;

  std::unique_ptr<SpatialIndex::IStorageManager> diskfile {
      SpatialIndex::StorageManager::createNewDiskStorageManager(filename,
                                                                4096)};

  // Configure RTree properties (see
  // https://libspatialindex.org/en/latest/overview.html#the-rtree-package)
  // Note: especially LeafCapacity and IndexCapacity are important for
  // performance
  /// @TODO: optimize these parameters based on the expected number of cells
  Tools::PropertySet ps;

  // Set FillFactor property
  Tools::Variant var;
  var.m_varType = Tools::VT_DOUBLE;
  var.m_val.dblVal = 0.7;
  ps.setProperty("FillFactor", var);

  // Set IndexCapacity property
  var.m_varType = Tools::VT_ULONG;
  var.m_val.ulVal = 16;
  ps.setProperty("IndexCapacity", var);

  // Set LeafCapacity property
  var.m_varType = Tools::VT_ULONG;
  var.m_val.ulVal = 16;
  ps.setProperty("LeafCapacity", var);

  // Set Dimension property (2D in our case)
  var.m_varType = Tools::VT_ULONG;
  var.m_val.ulVal = 2;
  ps.setProperty("Dimension", var);

  // Set TreeVariant property (R*-tree in our case)
  var.m_varType = Tools::VT_LONG;
  var.m_val.lVal = SpatialIndex::RTree::RV_RSTAR;
  ps.setProperty("TreeVariant", var);

  try {
    SpatialIndex::id_type indexId;
    std::unique_ptr<SpatialIndex::ISpatialIndex> tree {
        SpatialIndex::RTree::createAndBulkLoadNewRTree(
            SpatialIndex::RTree::BLM_STR, stream, *diskfile, ps, indexId)};

    // Clear temporary storage
    m_boxes.clear();
    m_boxes.shrink_to_fit();
  } catch (Tools::Exception& e) {
    throw std::runtime_error(std::string("Error building RTree: ") + e.what());
  }
}
