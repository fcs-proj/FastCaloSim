// Copyright (c) 2024 CERN for the benefit of the FastCaloSim project

#ifndef ISF_FASTCALOSIMEVENT_TFCSParametrizationAbsEtaSelectChain_h
#define ISF_FASTCALOSIMEVENT_TFCSParametrizationAbsEtaSelectChain_h

#include "FastCaloSim/Core/TFCSParametrizationEtaSelectChain.h"

class TFCSParametrizationAbsEtaSelectChain
    : public TFCSParametrizationEtaSelectChain
{
public:
  TFCSParametrizationAbsEtaSelectChain(const char* name = nullptr,
                                       const char* title = nullptr)
      : TFCSParametrizationEtaSelectChain(name, title) {};
  TFCSParametrizationAbsEtaSelectChain(
      const TFCSParametrizationAbsEtaSelectChain& ref)
      : TFCSParametrizationEtaSelectChain(ref) {};

  // selects on |extrapol->IDCaloBoundary_eta()|
  // return -1 if outside range
  virtual int get_bin(TFCSSimulationState&,
                      const TFCSTruthState* truth,
                      const TFCSExtrapolationState* extrapol) const override;
  virtual const std::string get_bin_text(int bin) const override;

private:
  ClassDefOverride(TFCSParametrizationAbsEtaSelectChain,
                   1)  // TFCSParametrizationAbsEtaSelectChain
};

#endif
