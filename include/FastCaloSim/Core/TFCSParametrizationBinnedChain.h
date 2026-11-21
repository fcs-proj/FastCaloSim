// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#pragma once

#include <cstddef>
#include <optional>
#include <string>

#include "FastCaloSim/Core/TFCSParametrizationChain.h"

class TFCSParametrizationBinnedChain : public TFCSParametrizationChain
{
public:
  explicit TFCSParametrizationBinnedChain(const char* name = nullptr,
                                          const char* title = nullptr)
      : TFCSParametrizationChain(name, title)
      , m_bin_start(1, 0)
  {
  }
  TFCSParametrizationBinnedChain(const TFCSParametrizationBinnedChain& ref) =
      default;
  TFCSParametrizationBinnedChain(TFCSParametrizationBinnedChain&& ref) =
      default;
  auto operator=(const TFCSParametrizationBinnedChain& ref)
      -> TFCSParametrizationBinnedChain& = default;
  auto operator=(TFCSParametrizationBinnedChain&& ref)
      -> TFCSParametrizationBinnedChain& = default;
  ~TFCSParametrizationBinnedChain() override = default;

  virtual void push_before_first_bin(TFCSParametrizationBase* param);
  virtual void push_back_in_bin(TFCSParametrizationBase* param,
                                std::size_t bin);

  virtual auto get_number_of_bins() const -> std::size_t
  {
    return m_bin_start.size() - 1;
  }

  /// this method should determine in derived classes which bin to simulate, so
  /// that the simulate method can call the appropriate TFCSParametrizationBase
  /// simulations; returns no value if no bin matches
  virtual auto get_bin(TFCSSimulationState& simulstate,
                       const TFCSTruthState* /*unused*/,
                       const TFCSExtrapolationState* /*unused*/) const
      -> std::optional<std::size_t>;
  virtual auto get_variable_text(TFCSSimulationState& /*unused*/,
                                 const TFCSTruthState* /*unused*/,
                                 const TFCSExtrapolationState* /*unused*/) const
      -> const std::string;
  /// print the range of a bin; for no bin or a bin out of range, print the
  /// allowed range
  virtual auto get_bin_text(std::optional<std::size_t> bin) const
      -> const std::string;

  auto simulate(TFCSSimulationState& simulstate,
                const TFCSTruthState* truth,
                const TFCSExtrapolationState* extrapol) const
      -> FCSReturnCode override;
  void Print(Option_t* option = "") const override;

protected:
  /// The bin to simulate if it is one of the bins of this chain
  auto find_bin(TFCSSimulationState& simulstate,
                const TFCSTruthState* truth,
                const TFCSExtrapolationState* extrapol) const
      -> std::optional<std::size_t>;
  /// The bin number, or "none"
  static auto bin_label(std::optional<std::size_t> bin) -> std::string;

  /// Contains the index where the TFCSParametrizationBase* instances to run for
  /// a given bin start. The last entry of the vector corresponds to the index
  /// from where on TFCSParametrizationBase* objects should be run again for all
  /// bins. This way one can loop over some instances for all bins, then only
  /// specific ones for one bin and at the end again over some for all bins
  std::vector<unsigned int> m_bin_start;

private:
  ClassDefOverride(TFCSParametrizationBinnedChain,
                   1)  // TFCSParametrizationBinnedChain
};
