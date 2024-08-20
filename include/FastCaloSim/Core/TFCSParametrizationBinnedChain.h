/*
  Copyright (C) 2002-2019 CERN for the benefit of the ATLAS collaboration
*/

#ifndef ISF_FASTCALOSIMEVENT_TFCSParametrizationBinnedChain_h
#define ISF_FASTCALOSIMEVENT_TFCSParametrizationBinnedChain_h

#include "FastCaloSim/Core/TFCSParametrizationChain.h"

class TFCSParametrizationBinnedChain : public TFCSParametrizationChain
{
public:
  TFCSParametrizationBinnedChain(const char* name = nullptr,
                                 const char* title = nullptr)
      : TFCSParametrizationChain(name, title)
      , m_bin_start(1, 0) {};
  TFCSParametrizationBinnedChain(const TFCSParametrizationBinnedChain& ref)
      : TFCSParametrizationChain(ref)
      , m_bin_start(ref.m_bin_start) {};

  virtual void push_before_first_bin(TFCSParametrizationBase* param);
  virtual void push_back_in_bin(TFCSParametrizationBase* param,
                                unsigned int bin);

  virtual unsigned int get_number_of_bins() const
  {
    return m_bin_start.size() - 1;
  };

  /// this method should determine in derived classes which bin to simulate, so
  /// that the simulate method can call the appropriate TFCSParametrizationBase
  /// simulations return -1 if no bin matches
  virtual int get_bin(TFCSSimulationState& simulstate,
                      const TFCSTruthState*,
                      const TFCSExtrapolationState*) const;
  virtual const std::string get_variable_text(
      TFCSSimulationState&,
      const TFCSTruthState*,
      const TFCSExtrapolationState*) const;
  /// print the range of a bin; for bin -1, print the allowed range
  virtual const std::string get_bin_text(int bin) const;

  virtual FCSReturnCode simulate(
      TFCSSimulationState& simulstate,
      const TFCSTruthState* truth,
      const TFCSExtrapolationState* extrapol) const override;
  void Print(Option_t* option = "") const override;

protected:
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

#endif
