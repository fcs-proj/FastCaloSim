// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#pragma once

#include "FastCaloSim/Core/TFCSParametrizationFloatSelectChain.h"
#include "FastCaloSim/Core/TFCSSimulationState.h"

class TFCSParametrizationEkinSelectChain
    : public TFCSParametrizationFloatSelectChain
{
public:
  TFCSParametrizationEkinSelectChain(const char* name = nullptr,
                                     const char* title = nullptr)
      : TFCSParametrizationFloatSelectChain(name, title)
  {
    reset_DoRandomInterpolation();
  };
  TFCSParametrizationEkinSelectChain(
      const TFCSParametrizationEkinSelectChain& ref)
      : TFCSParametrizationFloatSelectChain(ref)
  {
    reset_DoRandomInterpolation();
  };

  /// Status bit for Ekin Selection
  enum FCSEkinStatusBits
  {
    kDoRandomInterpolation =
        BIT(15)  ///< Set this bit in the TObject bit field if a random
                 ///< selection between neighboring Ekin bins should be done
  };

  auto DoRandomInterpolation() const -> bool
  {
    return TestBit(kDoRandomInterpolation);
  };
  void set_DoRandomInterpolation() { SetBit(kDoRandomInterpolation); };
  void reset_DoRandomInterpolation() { ResetBit(kDoRandomInterpolation); };

  using TFCSParametrizationFloatSelectChain::push_back_in_bin;
  virtual void push_back_in_bin(TFCSParametrizationBase* param);
  // selects on truth->Ekin()
  // return -1 if outside range
  auto get_bin(TFCSSimulationState&,
               const TFCSTruthState* truth,
               const TFCSExtrapolationState*) const -> int override;
  auto get_variable_text(TFCSSimulationState& simulstate,
                         const TFCSTruthState*,
                         const TFCSExtrapolationState*) const
      -> const std::string override;
  auto get_bin_text(int bin) const -> const std::string override;

protected:
  void recalc() override;

private:
  ClassDefOverride(TFCSParametrizationEkinSelectChain,
                   1)  // TFCSParametrizationEkinSelectChain
};
