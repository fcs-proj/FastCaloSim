// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#pragma once

#include "FastCaloSim/Core/TFCSParametrizationChain.h"

class TFCSParametrizationPDGIDSelectChain : public TFCSParametrizationChain
{
public:
  explicit TFCSParametrizationPDGIDSelectChain(const char* name = nullptr,
                                               const char* title = nullptr)
      : TFCSParametrizationChain(name, title)
  {
    reset_SimulateOnlyOnePDGID();
  }
  TFCSParametrizationPDGIDSelectChain(
      const TFCSParametrizationPDGIDSelectChain& ref) = default;
  TFCSParametrizationPDGIDSelectChain(
      TFCSParametrizationPDGIDSelectChain&& ref) = default;
  auto operator=(const TFCSParametrizationPDGIDSelectChain& ref)
      -> TFCSParametrizationPDGIDSelectChain& = default;
  auto operator=(TFCSParametrizationPDGIDSelectChain&& ref)
      -> TFCSParametrizationPDGIDSelectChain& = default;
  ~TFCSParametrizationPDGIDSelectChain() override = default;

  /// Status bit for PDGID Selection
  enum FCSPDGIDStatusBits
  {
    kSimulateOnlyOnePDGID = BIT(
        15)  ///< Set this bit in the TObject bit field if the PDGID selection
             ///< loop should be aborted after the first successful match
  };

  auto SimulateOnlyOnePDGID() const -> bool
  {
    return TestBit(kSimulateOnlyOnePDGID);
  }
  void set_SimulateOnlyOnePDGID() { SetBit(kSimulateOnlyOnePDGID); }
  void reset_SimulateOnlyOnePDGID() { ResetBit(kSimulateOnlyOnePDGID); }

  auto simulate(TFCSSimulationState& simulstate,
                const TFCSTruthState* truth,
                const TFCSExtrapolationState* extrapol) const
      -> FCSReturnCode override;

protected:
  void recalc() override;

private:
  ClassDefOverride(TFCSParametrizationPDGIDSelectChain,
                   1)  // TFCSParametrizationPDGIDSelectChain
};
