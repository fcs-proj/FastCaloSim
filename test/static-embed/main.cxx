// Copyright (c) 2025 CERN for the benefit of the FastCaloSim project

#include <cstdio>

// Provided by the shared library that absorbed the FastCaloSim archive.
auto fcs_count_missing_dictionaries() -> int;

auto main() -> int
{
  const int missing = fcs_count_missing_dictionaries();
  if (missing != 0) {
    std::printf("%d FastCaloSim dictionary(ies) missing after embedding\n",
                missing);
    return 1;
  }
  std::printf("all checked FastCaloSim dictionaries present after embedding\n");
  return 0;
}
