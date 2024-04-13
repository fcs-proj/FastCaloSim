/*
  Copyright (C) 2002-2019 CERN for the benefit of the ATLAS collaboration
*/

#ifndef ISF_FASTCALOSIMEVENT_TFCSEnergyRenormalization_h
#define ISF_FASTCALOSIMEVENT_TFCSEnergyRenormalization_h

#include "FastCaloSim/TFCSParametrization.h"

/** The class TFCSEnergyRenormalization ensures that the sum of cell energies in
   every calorimeter layer matches the output of energy simulation
*/

class TFCSEnergyRenormalization : public TFCSParametrization {
public:
  TFCSEnergyRenormalization(const char *name = nullptr,
                            const char *title = nullptr);
  virtual ~TFCSEnergyRenormalization();

  virtual bool is_match_Ekin_bin(int /*Ekin_bin*/) const override;
  virtual bool is_match_calosample(int /*calosample*/) const override;

  virtual FCSReturnCode
  simulate(TFCSSimulationState &simulstate, const TFCSTruthState * /*truth*/,
           const TFCSExtrapolationState * /*extrapol*/) const override;

private:
  ClassDefOverride(TFCSEnergyRenormalization, 1) // TFCSEnergyRenormalization
};

inline bool
TFCSEnergyRenormalization::is_match_Ekin_bin(int /*Ekin_bin*/) const {
  return true;
}

inline bool
TFCSEnergyRenormalization::is_match_calosample(int /*calosample*/) const {
  return true;
}

#endif
