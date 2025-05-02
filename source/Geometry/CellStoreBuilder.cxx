// Copyright (c) 2025 CERN for the benefit of the FastCaloSim project
#pragma once

#include <fstream>
#include <map>
#include <string>
#include <vector>

#include "FastCaloSim/Geometry/Cell.h"

/**
 * @brief A tool to write cells to disk for use with memory-mapped CellStore.
 *
 * Writes two files:
 * - <base_path>.data: flat binary file of serialized Cell objects
 * - <base_path>.index: sorted binary index (cell_id, offset)
 */
class CellStoreBuilder
{
public:
  void add_cell(const Cell& cell) { m_cells.push_back(cell); }

  void write(const std::string& base_path)
  {
    std::ofstream data_file(base_path + ".data", std::ios::binary);
    std::ofstream index_file(base_path + ".index", std::ios::binary);
    if (!data_file || !index_file) {
      throw std::runtime_error("Unable to open output files");
    }

    // Sort by ID for predictable access
    std::sort(m_cells.begin(),
              m_cells.end(),
              [](const Cell& a, const Cell& b) { return a.id() < b.id(); });

    for (const Cell& cell : m_cells) {
      uint64_t id = cell.id();
      uint64_t offset = static_cast<uint64_t>(data_file.tellp());
      index_file.write(reinterpret_cast<const char*>(&id), sizeof(uint64_t));
      index_file.write(reinterpret_cast<const char*>(&offset),
                       sizeof(uint64_t));
      data_file.write(reinterpret_cast<const char*>(&cell), sizeof(Cell));
    }
  }

private:
  std::vector<Cell> m_cells;
};
