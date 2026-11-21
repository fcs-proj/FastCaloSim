#include <print>

// -- Core FastCaloSim includes
#include "FastCaloSim/Core/TFCSParametrizationBase.h"

// -- Your own param gen includes
#include "FastCaloSimParam/YourClass.h"

auto main() -> int
{
  try {
    std::println("Running Parametrization Slice Creation...");
  } catch (...) {
    return 1;
  }
  return 0;
}
