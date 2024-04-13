/*
  Copyright (C) 2002-2018 CERN for the benefit of the ATLAS collaboration
*/

#ifndef ISF_FASTCALOSIMEVENT_TFCSParametrizationEtaSelectChain_h
#define ISF_FASTCALOSIMEVENT_TFCSParametrizationEtaSelectChain_h

#include "FastCaloSim/TFCSParametrizationFloatSelectChain.h"

class TFCSParametrizationEtaSelectChain
    : public TFCSParametrizationFloatSelectChain {
public:
  TFCSParametrizationEtaSelectChain(const char *name = nullptr,
                                    const char *title = nullptr)
      : TFCSParametrizationFloatSelectChain(name, title){};
  TFCSParametrizationEtaSelectChain(
      const TFCSParametrizationEtaSelectChain &ref)
      : TFCSParametrizationFloatSelectChain(ref){};

  using TFCSParametrizationFloatSelectChain::push_back_in_bin;
  virtual void push_back_in_bin(TFCSParametrizationBase *param);
  // selects on extrapol->IDCaloBoundary_eta()
  // return -1 if outside range
  virtual int get_bin(TFCSSimulationState &, const TFCSTruthState *truth,
                      const TFCSExtrapolationState *extrapol) const override;
  virtual const std::string
  get_variable_text(TFCSSimulationState &simulstate, const TFCSTruthState *,
                    const TFCSExtrapolationState *) const override;
  virtual const std::string get_bin_text(int bin) const override;

  static void unit_test(TFCSSimulationState *simulstate = nullptr,
                        TFCSTruthState *truth = nullptr,
                        TFCSExtrapolationState *extrapol = nullptr);

protected:
  virtual void recalc() override;

private:
  ClassDefOverride(TFCSParametrizationEtaSelectChain,
                   1) // TFCSParametrizationEtaSelectChain
};

#endif
