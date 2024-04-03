#include <string>

#include "FastCaloSim/FastCaloSim.hpp"

auto main() -> int
{
  auto const exported = exported_class {};

  return std::string("FastCaloSim") == exported.name() ? 0 : 1;
}
