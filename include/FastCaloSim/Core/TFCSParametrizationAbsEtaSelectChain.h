// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#pragma once

#include "FastCaloSim/Core/TFCSParametrizationEtaSelectChain.h"

class TFCSParametrizationAbsEtaSelectChain
    : public TFCSParametrizationEtaSelectChain
{
public:
  TFCSParametrizationAbsEtaSelectChain(const char* name = nullptr,
                                       const char* title = nullptr)
      : TFCSParametrizationEtaSelectChain(name, title)
  {
  }
  TFCSParametrizationAbsEtaSelectChain(
      const TFCSParametrizationAbsEtaSelectChain& ref) = default;

  // selects on |extrapol->IDCaloBoundary_eta()|
  // returns no value if outside range
  auto get_bin(TFCSSimulationState&,
               const TFCSTruthState* truth,
               const TFCSExtrapolationState* extrapol) const
      -> std::optional<std::size_t> override;
  auto get_bin_text(std::optional<std::size_t> bin) const
      -> const std::string override;

private:
  ClassDefOverride(TFCSParametrizationAbsEtaSelectChain,
                   1)  // TFCSParametrizationAbsEtaSelectChain
};
