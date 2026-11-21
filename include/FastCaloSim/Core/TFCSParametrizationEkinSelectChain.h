// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#pragma once

#include "FastCaloSim/Core/TFCSParametrizationFloatSelectChain.h"
#include "FastCaloSim/Core/TFCSSimulationState.h"

// Custom copy constructor, moves fall back to it
// NOLINTNEXTLINE(cppcoreguidelines-special-member-functions,hicpp-special-member-functions)
class TFCSParametrizationEkinSelectChain
    : public TFCSParametrizationFloatSelectChain
{
public:
  explicit TFCSParametrizationEkinSelectChain(const char* name = nullptr,
                                              const char* title = nullptr)
      : TFCSParametrizationFloatSelectChain(name, title)
  {
    reset_DoRandomInterpolation();
  }
  TFCSParametrizationEkinSelectChain(
      const TFCSParametrizationEkinSelectChain& ref)
      : TFCSParametrizationFloatSelectChain(ref)
  {
    reset_DoRandomInterpolation();
  }

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
  }
  void set_DoRandomInterpolation() { SetBit(kDoRandomInterpolation); }
  void reset_DoRandomInterpolation() { ResetBit(kDoRandomInterpolation); }

  using TFCSParametrizationFloatSelectChain::push_back_in_bin;
  virtual void push_back_in_bin(TFCSParametrizationBase* param);
  // selects on truth->Ekin()
  // returns no value if outside range
  auto get_bin(TFCSSimulationState& /*simulstate*/,
               const TFCSTruthState* truth,
               const TFCSExtrapolationState* /*unused*/) const
      -> std::optional<std::size_t> override;
  auto get_variable_text(TFCSSimulationState& simulstate,
                         const TFCSTruthState* /*truth*/,
                         const TFCSExtrapolationState* /*unused*/) const
      -> const std::string override;
  auto get_bin_text(std::optional<std::size_t> bin) const
      -> const std::string override;

protected:
  void recalc() override;

private:
  ClassDefOverride(TFCSParametrizationEkinSelectChain,
                   1)  // TFCSParametrizationEkinSelectChain
};
