// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#ifndef ISF_FASTCALOSIMEVENT_TFCSEnergyRenormalization_h
#define ISF_FASTCALOSIMEVENT_TFCSEnergyRenormalization_h

#include "FastCaloSim/Core/TFCSParametrization.h"

class CaloGeo;

/** The class TFCSEnergyRenormalization ensures that the sum of cell energies in
   every calorimeter layer matches the output of energy simulation
*/

class TFCSEnergyRenormalization : public TFCSParametrization
{
public:
  TFCSEnergyRenormalization(const char* name = nullptr,
                            const char* title = nullptr,
                            CaloGeo* geo = nullptr);

  void set_geometry(CaloGeo* geo) override { m_geo = geo; };
  ~TFCSEnergyRenormalization() override;

  auto is_match_Ekin_bin(int /*Ekin_bin*/) const -> bool override;
  auto is_match_calosample(int /*calosample*/) const -> bool override;

  auto simulate(TFCSSimulationState& simulstate,
                const TFCSTruthState* /*truth*/,
                const TFCSExtrapolationState* /*extrapol*/) const
      -> FCSReturnCode override;

protected:
  CaloGeo* m_geo;  //! do not persistify

private:
  ClassDefOverride(TFCSEnergyRenormalization, 1)  // TFCSEnergyRenormalization
};

inline auto TFCSEnergyRenormalization::is_match_Ekin_bin(int /*Ekin_bin*/) const
    -> bool
{
  return true;
}

inline auto TFCSEnergyRenormalization::is_match_calosample(
    int /*calosample*/) const -> bool
{
  return true;
}

#endif
