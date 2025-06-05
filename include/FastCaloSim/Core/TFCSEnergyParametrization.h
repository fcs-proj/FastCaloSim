// Copyright (c) 2024 CERN for the benefit of the FastCaloSim project

#ifndef ISF_FASTCALOSIMEVENT_TFCSEnergyParametrization_h
#define ISF_FASTCALOSIMEVENT_TFCSEnergyParametrization_h

#include "FastCaloSim/Core/TFCSParametrization.h"

namespace FastCaloSim::Core
{
class TFCSEnergyParametrization : public TFCSParametrization
{
public:
  TFCSEnergyParametrization(const char* name = nullptr,
                            const char* title = nullptr);

  virtual bool is_match_Ekin_bin(int /*Ekin_bin*/) const override
  {
    return true;
  };
  virtual bool is_match_calosample(int /*calosample*/) const override
  {
    return true;
  };

  // return number of energy parametrization bins
  virtual int n_bins() const { return 0; };

private:
  ClassDefOverride(TFCSEnergyParametrization, 1)  // TFCSEnergyParametrization
};
}  // namespace FastCaloSim::Core

#endif
