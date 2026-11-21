// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#pragma once

#include <concepts>

namespace FastCaloSim
{
/// Exact comparison of two floating-point numbers, for the cases where exact
/// equality is meant: a value that is exactly zero (e.g. set to zero to mark a
/// special case) or ties between values computed in the same way.
/// -Wfloat-equal is suppressed only here, so that it still reports unintended
/// exact comparisons everywhere else.
template<std::floating_point T>
constexpr auto exactly_equal(T a, T b) -> bool
{
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-equal"
  return a == b;
#pragma GCC diagnostic pop
}
}  // namespace FastCaloSim
