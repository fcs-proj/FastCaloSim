// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#pragma once

#include <algorithm>
#include <bit>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <span>
#include <type_traits>
#include <vector>

#include "FastCaloSim/Definitions/Index.h"

namespace FastCaloSim
{
/// Guide table (Chen and Asau, 1974) for sampling from a cumulative
/// distribution c[0..n-1] (non-decreasing) with a uniform random number u.
///
/// upper_bound(c, u) returns exactly what std::ranges::upper_bound returns
/// (the index of the first element with u < c[i], or n), in O(1) on average
/// instead of O(log n): the range of u is divided into 2^b equal buckets, at
/// least as many as there are bins, and the table stores the result for the
/// lower edge of every bucket. The result for a u inside a bucket lies between
/// the results for the bucket's two edges, so only that part of c is searched,
/// on average less than one element. The bucket of u and its edges are powers
/// of two apart and computed without rounding, so the result is the same as
/// that of the full search for every u.
///
/// Two kinds of cumulative distributions are supported:
///   - float values in [0, 1], searched with a double u in [0, 1); for any
///     other u (negative, >= 1, NaN) the full search is used
///   - 32-bit unsigned integers, searched with an unsigned 32-bit u
/// The table is built only for sorted c without NaN; otherwise, and when
/// upper_bound() gets a c of another size than the table was built for, the
/// full search is used, so the result is always that of the full search.
///
/// build() must be called whenever c changes; upper_bound() only reads the
/// table and can be called from several threads.
template<typename T>
  requires std::same_as<T, float> || std::same_as<T, std::uint32_t>
class GuideTable
{
public:
  /// The type of the random number: floats are compared in double precision
  using rnd_type = std::conditional_t<std::same_as<T, float>, double, T>;

  void build(std::span<const T> cumulative)
  {
    m_first.clear();
    m_first.shrink_to_fit();
    m_bits = 0;
    m_nbuckets = 1;
    m_size = cumulative.size();
    if (cumulative.empty()
        || cumulative.size() > std::numeric_limits<std::uint32_t>::max())
    {
      return;
    }
    if constexpr (std::same_as<T, float>) {
      if (std::ranges::any_of(cumulative,
                              [](T c) -> bool { return std::isnan(c); }))
      {
        return;
      }
    }
    if (!std::ranges::is_sorted(cumulative)) {
      return;
    }

    // At least as many buckets as bins, at most 2^24
    m_bits = std::min(static_cast<unsigned>(std::bit_width(m_size - 1)), 24U);
    const std::size_t nbuckets = std::size_t {1} << m_bits;
    m_nbuckets = static_cast<double>(nbuckets);
    m_first.resize(nbuckets + 1);
    // One pass over the bucket edges and c together
    std::size_t i = 0;
    for (std::size_t k = 0; k <= nbuckets; ++k) {
      while (i < m_size && !before(k, cumulative[i])) {
        ++i;
      }
      m_first[k] = static_cast<std::uint32_t>(i);
    }
  }

  /// Same result as std::ranges::upper_bound(cumulative, u) (floats compared
  /// in double precision), given the cumulative distribution the table was
  /// built for
  auto upper_bound(std::span<const T> cumulative, rnd_type u) const
      -> std::size_t
  {
    if constexpr (std::same_as<T, float>) {
      if (m_first.empty() || cumulative.size() != m_size
          || !(u >= 0.0 && u < 1.0))
      {
        return full_search(cumulative, u);
      }
      // u * 2^b is exact, the bucket is its integer part
      const auto k = static_cast<std::size_t>(u * m_nbuckets);
      return search(cumulative, k, u);
    } else {
      if (m_first.empty() || cumulative.size() != m_size) {
        return full_search(cumulative, u);
      }
      const auto k =
          static_cast<std::size_t>((std::uint64_t {u} << m_bits) >> 32U);
      return search(cumulative, k, u);
    }
  }

  auto empty() const -> bool { return m_first.empty(); }

  /// Number of table entries (buckets + 1)
  auto size() const -> std::size_t { return m_first.size(); }

private:
  static auto full_search(std::span<const T> cumulative, rnd_type u)
      -> std::size_t
  {
    const auto it = std::ranges::upper_bound(
        cumulative,
        u,
        std::ranges::less {},
        [](T c) -> rnd_type { return static_cast<rnd_type>(c); });
    return index_of(cumulative, it);
  }

  /// The result lies between the results for the two edges of bucket k
  auto search(std::span<const T> cumulative, std::size_t k, rnd_type u) const
      -> std::size_t
  {
    const auto it = std::ranges::upper_bound(
        cumulative.begin() + m_first[k],
        cumulative.begin() + m_first[k + 1],
        u,
        std::ranges::less {},
        [](T c) -> rnd_type { return static_cast<rnd_type>(c); });
    return index_of(cumulative, it);
  }

  /// Is the lower edge of bucket k smaller than the value c?
  auto before(std::size_t k, T c) const -> bool
  {
    if constexpr (std::same_as<T, float>) {
      // k / 2^b is exact
      return static_cast<double>(k) / m_nbuckets < static_cast<double>(c);
    } else {
      return (std::uint64_t {k} << (32U - m_bits)) < std::uint64_t {c};
    }
  }

  /// m_first[k]: result of the search for the lower edge of bucket k
  std::vector<std::uint32_t> m_first;
  unsigned m_bits = 0;
  /// 2^m_bits, the factor between u in [0, 1) and its bucket (floats)
  double m_nbuckets = 1;
  std::size_t m_size = 0;
};
}  // namespace FastCaloSim
