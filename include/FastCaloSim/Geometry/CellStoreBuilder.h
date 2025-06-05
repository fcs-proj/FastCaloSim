#include <algorithm>
#include <fstream>
#include <stdexcept>
#include <vector>

#include "FastCaloSim/Geometry/Cell.h"

namespace FastCaloSim::Geometry
{

class CellStoreBuilder
{
public:
  void add_cell(const Cell& cell) { m_cells.push_back(cell.raw()); }

  void write(const std::string& base_path)
  {
    std::ofstream data_file(base_path + ".data", std::ios::binary);
    std::ofstream index_file(base_path + ".index", std::ios::binary);
    if (!data_file || !index_file) {
      throw std::runtime_error("Unable to open output files");
    }

    // Sort by ID, needed for binary search in lookup
    std::sort(m_cells.begin(),
              m_cells.end(),
              [](const CellData& a, const CellData& b)
              { return a.m_id < b.m_id; });

    for (const CellData& data : m_cells) {
      uint64_t id = data.m_id;

      std::streampos pos = data_file.tellp();
      if (pos == std::streampos(-1)) {
        throw std::runtime_error("tellp() failed");
      }

      uint64_t offset = static_cast<uint64_t>(pos);
      index_file.write(reinterpret_cast<const char*>(&id), sizeof(uint64_t));
      index_file.write(reinterpret_cast<const char*>(&offset),
                       sizeof(uint64_t));
      data_file.write(reinterpret_cast<const char*>(&data), sizeof(CellData));
    }
    // Free memory after writing
    m_cells.clear();
    m_cells.shrink_to_fit();
    data_file.close();
    index_file.close();
  }

private:
  std::vector<CellData> m_cells;
};
}  // namespace FastCaloSim::Geometry
