// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#pragma once

#include <cstddef>
#include <ranges>

namespace FastCaloSim
{
/// The index of an iterator into a random access range, e.g. the result of
/// std::upper_bound on a std::vector
template<std::ranges::random_access_range Range>
constexpr auto index_of(const Range& range,
                        std::ranges::iterator_t<const Range> it) -> std::size_t
{
  return static_cast<std::size_t>(it - std::ranges::begin(range));
}
}  // namespace FastCaloSim
