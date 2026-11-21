// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#pragma once

#include "FastCaloSim/Core/TFCSParametrizationFloatSelectChain.h"

class TFCSParametrizationEtaSelectChain
    : public TFCSParametrizationFloatSelectChain
{
public:
  TFCSParametrizationEtaSelectChain(const char* name = nullptr,
                                    const char* title = nullptr)
      : TFCSParametrizationFloatSelectChain(name, title)
  {
  }
  TFCSParametrizationEtaSelectChain(
      const TFCSParametrizationEtaSelectChain& ref) = default;

  using TFCSParametrizationFloatSelectChain::push_back_in_bin;
  virtual void push_back_in_bin(TFCSParametrizationBase* param);
  // selects on extrapol->IDCaloBoundary_eta()
  // returns no value if outside range
  auto get_bin(TFCSSimulationState&,
               const TFCSTruthState* truth,
               const TFCSExtrapolationState* extrapol) const
      -> std::optional<std::size_t> override;
  auto get_variable_text(TFCSSimulationState& simulstate,
                         const TFCSTruthState*,
                         const TFCSExtrapolationState*) const
      -> const std::string override;
  auto get_bin_text(std::optional<std::size_t> bin) const
      -> const std::string override;

protected:
  void recalc() override;

private:
  ClassDefOverride(TFCSParametrizationEtaSelectChain,
                   1)  // TFCSParametrizationEtaSelectChain
};
