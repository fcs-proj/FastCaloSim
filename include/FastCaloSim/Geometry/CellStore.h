// Copyright (c) 2025 CERN for the benefit of the FastCaloSim project
#pragma once
#include <cstdint>
#include <cstring>
#include <fstream>
#include <memory>
#include <stdexcept>
#include <string>

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
  /// @brief Constructor
  CellStore() = default;

  ~CellStore() { unmap(); }

  // CellStore owns mmap'd regions, so it must not be copied (a shallow copy
  // would double-munmap). It is movable: ownership of the mappings transfers
  // and the source is left empty.
  CellStore(const CellStore&) = delete;
  auto operator=(const CellStore&) -> CellStore& = delete;

  CellStore(CellStore&& other) noexcept { *this = std::move(other); }

  auto operator=(CellStore&& other) noexcept -> CellStore&
  {
    if (this != &other) {
      unmap();
      m_data = other.m_data;
      m_file_size = other.m_file_size;
      m_n_cells = other.m_n_cells;
      m_index_data = other.m_index_data;
      m_index_size = other.m_index_size;
      other.m_data = nullptr;
      other.m_index_data = nullptr;
      other.m_file_size = 0;
      other.m_index_size = 0;
      other.m_n_cells = 0;
    }
    return *this;
  }

  /// @brief Load cell data and index from the specified base path
  /// @param base_path The base path for the data and index files
  /// @param cache_size_bytes Unused: the OS page cache supersedes the former
  ///        user-space LRU cache. Kept for API compatibility.
  void load(const std::string& base_path,
            [[maybe_unused]] size_t cache_size_bytes = 10 * 1024 * 1024)
  {
    // Release any mappings from a previous load() so reloading does not leak
    unmap();

    std::string data_path = base_path + ".data";
    std::string index_path = base_path + ".index";

    // Verify the index file exists/opens before mapping anything
    std::ifstream index_file(index_path, std::ios::binary);
    if (!index_file) {
      throw std::runtime_error("Failed to open index file: " + index_path);
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

    // Memory map the index file for binary search
    fd = ::open(index_path.c_str(), O_RDONLY);
    if (fd == -1) {
      munmap(m_data, m_file_size);
      m_data = nullptr;
      throw std::runtime_error("Failed to open index file: " + index_path);
    }

    if (fstat(fd, &st) == -1) {
      ::close(fd);
      munmap(m_data, m_file_size);
      m_data = nullptr;
      throw std::runtime_error("Failed to stat index file: " + index_path);
    }

    m_index_size = st.st_size;
    m_index_data = mmap(nullptr, m_index_size, PROT_READ, MAP_PRIVATE, fd, 0);
    ::close(fd);

    if (m_index_data == MAP_FAILED) {
      munmap(m_data, m_file_size);
      m_data = nullptr;
      throw std::runtime_error("Failed to mmap index file: " + index_path);
    }
  }

  /// @brief Return the cell with a given ID.
  /// @note The returned reference aliases a per-thread scratch Cell that is
  ///       overwritten by the next get()/get_at_index() call on the same
  ///       thread. Copy the Cell if you need to hold more than one at a time.
  auto get(uint64_t id) const -> const Cell&
  {
    if (m_data == nullptr || m_index_data == nullptr) {
      throw std::runtime_error("CellStore not loaded - call load() first");
    }
    return getCellAtOffset(findOffsetById(id));
  }

  /// @brief Return the cell at a given dense index.
  /// @note The returned reference aliases a per-thread scratch Cell that is
  ///       overwritten by the next get()/get_at_index() call on the same
  ///       thread. Copy the Cell if you need to hold more than one at a time.
  auto get_at_index(size_t idx) const -> const Cell&
  {
    if (!m_data) {
      throw std::runtime_error("CellStore not loaded - call load() first");
    }

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
  /// @brief Find the offset for a given cell ID using binary search in mmap'd
  /// index
  auto findOffsetById(uint64_t id) const -> uint64_t
  {
    // Cast the memory-mapped index to an array of id-offset pairs
    struct IdOffsetPair
    {
      uint64_t id;
      uint64_t offset;
    };

    const IdOffsetPair* pairs = static_cast<const IdOffsetPair*>(m_index_data);
    size_t num_pairs = m_index_size / sizeof(IdOffsetPair);

    // Binary search (we know the index is sorted by ID thanks to
    // CellStoreBuilder)
    size_t left = 0;
    size_t right = num_pairs;

    while (left < right) {
      size_t mid = left + (right - left) / 2;
      if (pairs[mid].id == id) {
        return pairs[mid].offset;
      } else if (pairs[mid].id < id) {
        left = mid + 1;
      } else {
        right = mid;
      }
    }

    throw std::runtime_error("Cell ID not found in index: "
                             + std::to_string(id));
  }

  /// @brief Release any memory-mapped regions
  ///
  /// Safe to call repeatedly: pointers and sizes are reset so that a
  /// subsequent call (or a destructor after a move) is a no-op.
  void unmap() noexcept
  {
    if (m_data != nullptr) {
      munmap(m_data, m_file_size);
      m_data = nullptr;
    }
    if (m_index_data != nullptr) {
      munmap(m_index_data, m_index_size);
      m_index_data = nullptr;
    }
    m_file_size = 0;
    m_index_size = 0;
    m_n_cells = 0;
  }

  /// @brief Get a cell at the specified offset in the data file
  auto getCellAtOffset(uint64_t offset) const -> const Cell&
  {
    const CellData* cell_data = reinterpret_cast<const CellData*>(
        static_cast<const char*>(m_data) + offset);

    // Use thread-local Cell to avoid lifetime issues
    thread_local Cell tmp_cell;
    std::memcpy(&tmp_cell.raw(), cell_data, sizeof(CellData));
    return tmp_cell;
  }

private:
  /// @brief Pointer to the memory-mapped data file containing serialized
  /// CellData
  void* m_data = nullptr;
  /// @brief Size of the memory-mapped file in bytes
  size_t m_file_size = 0;
  /// @brief Total number of cells stored (used for iteration)
  size_t m_n_cells = 0;

  /// @brief Pointer to memory-mapped index file
  void* m_index_data = nullptr;
  /// @brief Size of the index file in bytes
  size_t m_index_size = 0;
};
