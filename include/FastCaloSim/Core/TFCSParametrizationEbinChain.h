// Copyright (c) 2024 CERN for the benefit of the FastCaloSim project

#ifndef TFCSParametrizationEbinChain_h
#define TFCSParametrizationEbinChain_h

#include "FastCaloSim/Core/TFCSParametrizationBinnedChain.h"
#include "FastCaloSim/Core/TFCSSimulationState.h"

namespace FastCaloSim::Core
{
class TFCSParametrizationEbinChain : public TFCSParametrizationBinnedChain
{
public:
  TFCSParametrizationEbinChain(const char* name = nullptr,
                               const char* title = nullptr)
      : TFCSParametrizationBinnedChain(name, title) {};
  TFCSParametrizationEbinChain(const TFCSParametrizationEbinChain& ref)
      : TFCSParametrizationBinnedChain(ref) {};

  /// current convention is to start Ebin counting at 1, to be updated to start
  /// counting with 0
  virtual int get_bin(TFCSSimulationState& simulstate,
                      const TFCSTruthState*,
                      const TFCSExtrapolationState*) const override
  {
    return simulstate.Ebin();
  };
  virtual const std::string get_variable_text(
      TFCSSimulationState& simulstate,
      const TFCSTruthState*,
      const TFCSExtrapolationState*) const override;

private:
  ClassDefOverride(TFCSParametrizationEbinChain,
                   1)  // TFCSParametrizationEbinChain
};
}  // namespace FastCaloSim::Core

#endif
