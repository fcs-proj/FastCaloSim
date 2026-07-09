// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#ifndef ISF_FASTCALOSIMEVENT_TFCSParametrizationEbinChain_h
#define ISF_FASTCALOSIMEVENT_TFCSParametrizationEbinChain_h

#include "FastCaloSim/Core/TFCSParametrizationBinnedChain.h"
#include "FastCaloSim/Core/TFCSSimulationState.h"

class TFCSParametrizationEbinChain : public TFCSParametrizationBinnedChain
{
public:
  TFCSParametrizationEbinChain(const char* name = nullptr,
                               const char* title = nullptr)
      : TFCSParametrizationBinnedChain(name, title) {};
  TFCSParametrizationEbinChain(const TFCSParametrizationEbinChain& ref) =
      default;

  /// current convention is to start Ebin counting at 1, to be updated to start
  /// counting with 0
  auto get_bin(TFCSSimulationState& simulstate,
               const TFCSTruthState*,
               const TFCSExtrapolationState*) const -> int override
  {
    return simulstate.Ebin();
  };
  auto get_variable_text(TFCSSimulationState& simulstate,
                         const TFCSTruthState*,
                         const TFCSExtrapolationState*) const
      -> const std::string override;

private:
  ClassDefOverride(TFCSParametrizationEbinChain,
                   1)  // TFCSParametrizationEbinChain
};

#endif
