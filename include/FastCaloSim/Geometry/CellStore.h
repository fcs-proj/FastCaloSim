// Copyright (c) 2025 CERN for the benefit of the FastCaloSim project

#pragma once

#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

#include "FastCaloSim/Geometry/Cell.h"

/**
 * @brief A memory-mapped read-only store of CellData, indexed by ID.
 */
class CellStore
{
public:
  CellStore() = default;

  ~CellStore()
  {
    if (m_data != nullptr) {
      munmap(m_data, m_file_size);
    }
  }

  void load(const std::string& base_path)
  {
    std::string data_path = base_path + ".data";
    std::string index_path = base_path + ".index";

    // Load index file
    std::ifstream index_file(index_path, std::ios::binary);
    if (!index_file) {
      throw std::runtime_error("Failed to open index file: " + index_path);
    }

    for (;;) {
      uint64_t id, offset;
      if (!index_file.read(reinterpret_cast<char*>(&id), sizeof(id)))
        break;
      if (!index_file.read(reinterpret_cast<char*>(&offset), sizeof(offset))) {
        throw std::runtime_error("Corrupted index file (offset missing)");
      }
      m_index[id] = offset;
    }

    // Open and memory map the data file
    int fd = ::open(data_path.c_str(), O_RDONLY);
    if (fd == -1) {
      throw std::runtime_error("Failed to open data file: " + data_path);
    }

    struct stat st;
    if (fstat(fd, &st) == -1) {
      ::close(fd);
      throw std::runtime_error("Failed to stat data file: " + data_path);
    }

    m_file_size = st.st_size;
    m_data = mmap(nullptr, m_file_size, PROT_READ, MAP_PRIVATE, fd, 0);
    ::close(fd);

    if (m_data == MAP_FAILED) {
      throw std::runtime_error("Failed to mmap data file: " + data_path);
    }

    m_n_cells = m_file_size / sizeof(CellData);
  }

  auto get(uint64_t id) const -> const Cell&
  {
    auto it = m_index.find(id);
    if (it == m_index.end()) {
      throw std::runtime_error("Cell ID not found in index: "
                               + std::to_string(id));
      return m_invalid_cell;
    }

    const CellData* cell_data = reinterpret_cast<const CellData*>(
        static_cast<const char*>(m_data) + it->second);
    return *reinterpret_cast<const Cell*>(cell_data);
  }

  auto get_at_index(size_t idx) const -> const Cell&
  {
    if (idx >= m_n_cells) {
      throw std::runtime_error("Index out of bounds: " + std::to_string(idx));
    }

    // Calculate raw pointer to the bytes of CellData
    const char* cell_bytes =
        static_cast<const char*>(m_data) + idx * sizeof(CellData);

    // Use thread-local Cell to avoid lifetime issues
    thread_local Cell tmp_cell;
    std::memcpy(&tmp_cell.raw(), cell_bytes, sizeof(CellData));

    return tmp_cell;
  }

  auto size() const -> size_t { return m_n_cells; }

private:
  /// @brief Pointer to the memory-mapped data file containing serialized
  /// CellData
  void* m_data = nullptr;
  /// @brief Size of the memory-mapped file in bytes
  size_t m_file_size = 0;
  /// @brief Total number of cells stored (used for iteration)
  size_t m_n_cells = 0;
  /// @brief Mapping from cell ID to byte offset in the data file
  std::unordered_map<uint64_t, uint64_t> m_index;
  /// @brief Static fallback Cell instance returned for invalid or missing IDs
  static inline Cell m_invalid_cell {};
};
