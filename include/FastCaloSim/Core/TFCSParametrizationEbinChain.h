/*
  Copyright (C) 2002-2018 CERN for the benefit of the ATLAS collaboration
*/

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

  static void unit_test(TFCSSimulationState* simulstate = nullptr,
                        const TFCSTruthState* truth = nullptr,
                        const TFCSExtrapolationState* extrapol = nullptr);

private:
  ClassDefOverride(TFCSParametrizationEbinChain,
                   1)  // TFCSParametrizationEbinChain
};

#endif
