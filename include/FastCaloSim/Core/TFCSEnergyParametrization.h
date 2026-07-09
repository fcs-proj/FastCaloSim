// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#ifndef ISF_FASTCALOSIMEVENT_TFCSEnergyParametrization_h
#define ISF_FASTCALOSIMEVENT_TFCSEnergyParametrization_h

#include "FastCaloSim/Core/TFCSParametrization.h"

class TFCSEnergyParametrization : public TFCSParametrization
{
public:
  TFCSEnergyParametrization(const char* name = nullptr,
                            const char* title = nullptr);

  auto is_match_Ekin_bin(int /*Ekin_bin*/) const -> bool override
  {
    return true;
  };
  auto is_match_calosample(int /*calosample*/) const -> bool override
  {
    return true;
  };

  // return number of energy parametrization bins
  virtual auto n_bins() const -> int { return 0; };

private:
  ClassDefOverride(TFCSEnergyParametrization, 1)  // TFCSEnergyParametrization
};

#endif
