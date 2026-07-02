// Copyright (c) 2025 CERN for the benefit of the FastCaloSim project

#include <cstdio>

#include "TClass.h"

// A representative spread of FastCaloSim classes that are streamed back from
// the parametrization file at runtime. Their ROOT dictionaries live in
// FastCaloSim's generated dictionary translation unit; if that unit were
// dropped when the static archive was absorbed into this shared library, the
// lookups below would return null.
namespace
{
const char* const kClasses[] = {
    "TFCSParametrizationBase",
    "TFCSParametrization",
    "TFCSInvisibleParametrization",
    "TFCSEnergyInterpolationLinear",
    "TFCSPCAEnergyParametrization",
    "TFCSParametrizationEbinChain",
    "TFCSTruthState",
    "TFCSExtrapolationState",
    "TFCSSimulationState",
};
}  // namespace

// Returns the number of expected dictionaries that are missing (0 == success).
auto fcs_count_missing_dictionaries() -> int
{
  int missing = 0;
  for (const char* name : kClasses) {
    if (TClass::GetClass(name) == nullptr) {
      std::printf("missing dictionary: %s\n", name);
      ++missing;
    }
  }
  return missing;
}
