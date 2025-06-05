// Copyright (c) 2024 CERN for the benefit of the FastCaloSim project

#ifndef TFCSEnergyRenormalization_h
#define TFCSEnergyRenormalization_h

#include "FastCaloSim/Core/TFCSParametrization.h"

namespace FastCaloSim::Geometry
{
class CaloGeo;
}
/** The class TFCSEnergyRenormalization ensures that the sum of cell energies in
   every calorimeter layer matches the output of energy simulation
*/
namespace FastCaloSim::Core
{
class TFCSEnergyRenormalization : public TFCSParametrization
{
public:
  TFCSEnergyRenormalization(const char* name = nullptr,
                            const char* title = nullptr,
                            CaloGeo* geo = nullptr);

  virtual void set_geometry(CaloGeo* geo) override { m_geo = geo; };
  virtual ~TFCSEnergyRenormalization();

  virtual bool is_match_Ekin_bin(int /*Ekin_bin*/) const override;
  virtual bool is_match_calosample(int /*calosample*/) const override;

  virtual FCSReturnCode simulate(
      TFCSSimulationState& simulstate,
      const TFCSTruthState* /*truth*/,
      const TFCSExtrapolationState* /*extrapol*/) const override;

protected:
  CaloGeo* m_geo;  //! do not persistify

private:
  ClassDefOverride(TFCSEnergyRenormalization, 1)  // TFCSEnergyRenormalization
};

inline bool TFCSEnergyRenormalization::is_match_Ekin_bin(int /*Ekin_bin*/) const
{
  return true;
}

inline bool TFCSEnergyRenormalization::is_match_calosample(
    int /*calosample*/) const
{
  return true;
}
}  // namespace FastCaloSim::Core

#endif
