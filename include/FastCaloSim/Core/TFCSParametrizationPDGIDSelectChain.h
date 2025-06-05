// Copyright (c) 2024 CERN for the benefit of the FastCaloSim project

#ifndef TFCSParametrizationPDGIDSelectChain_h
#define TFCSParametrizationPDGIDSelectChain_h

#include "FastCaloSim/Core/TFCSParametrizationChain.h"

namespace FastCaloSim::Core
{

class TFCSParametrizationPDGIDSelectChain : public TFCSParametrizationChain
{
public:
  TFCSParametrizationPDGIDSelectChain(const char* name = nullptr,
                                      const char* title = nullptr)
      : TFCSParametrizationChain(name, title)
  {
    reset_SimulateOnlyOnePDGID();
  };
  TFCSParametrizationPDGIDSelectChain(
      const TFCSParametrizationPDGIDSelectChain& ref)
      : TFCSParametrizationChain(ref)
  {
    reset_SimulateOnlyOnePDGID();
  };

  /// Status bit for PDGID Selection
  enum FCSPDGIDStatusBits
  {
    kSimulateOnlyOnePDGID = BIT(
        15)  ///< Set this bit in the TObject bit field if the PDGID selection
             ///< loop should be aborted after the first successful match
  };

  bool SimulateOnlyOnePDGID() const { return TestBit(kSimulateOnlyOnePDGID); };
  void set_SimulateOnlyOnePDGID() { SetBit(kSimulateOnlyOnePDGID); };
  void reset_SimulateOnlyOnePDGID() { ResetBit(kSimulateOnlyOnePDGID); };

  virtual FCSReturnCode simulate(
      TFCSSimulationState& simulstate,
      const TFCSTruthState* truth,
      const TFCSExtrapolationState* extrapol) const override;

protected:
  virtual void recalc() override;

private:
  ClassDefOverride(TFCSParametrizationPDGIDSelectChain,
                   1)  // TFCSParametrizationPDGIDSelectChain
};
}  // namespace FastCaloSim::Core
#endif
