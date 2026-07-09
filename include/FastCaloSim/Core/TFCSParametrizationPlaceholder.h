// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#ifndef ISF_FASTCALOSIMEVENT_TFCSParametrizationPlaceholder_h
#define ISF_FASTCALOSIMEVENT_TFCSParametrizationPlaceholder_h

#include "FastCaloSim/Core/TFCSParametrizationBase.h"

class TFCSParametrizationPlaceholder : public TFCSParametrizationBase
{
public:
  TFCSParametrizationPlaceholder(const char* name = nullptr,
                                 const char* title = nullptr)
      : TFCSParametrizationBase(name, title) {};

  auto is_match_Ekin_bin(int /*Ekin_bin*/) const -> bool override
  {
    return true;
  };
  auto is_match_calosample(int /*calosample*/) const -> bool override
  {
    return true;
  };

  auto simulate(TFCSSimulationState& simulstate,
                const TFCSTruthState* truth,
                const TFCSExtrapolationState* extrapol) const
      -> FCSReturnCode override;

private:
  ClassDefOverride(TFCSParametrizationPlaceholder,
                   1)  // TFCSParametrizationPlaceholder
};

#endif
