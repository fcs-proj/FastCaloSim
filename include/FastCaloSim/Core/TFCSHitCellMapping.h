// Copyright (c) 2024 CERN for the benefit of the FastCaloSim project

#ifndef TFCSHitCellMapping_h
#define TFCSHitCellMapping_h

#include "FastCaloSim/Core/TFCSLateralShapeParametrizationHitBase.h"
#include "FastCaloSim/Geometry/Cell.h"

namespace FastCaloSim::Geometry
{
class CaloGeo;
}
namespace FastCaloSim::Core
{
class TFCSHitCellMapping : public TFCSLateralShapeParametrizationHitBase
{
public:
  TFCSHitCellMapping(const char* name = nullptr,
                     const char* title = nullptr,
                     CaloGeo* geo = nullptr);

  virtual void set_geometry(CaloGeo* geo) override { m_geo = geo; };
  CaloGeo* get_geometry() { return m_geo; };

  /// fills all hits into calorimeter cells
  virtual FCSReturnCode simulate_hit(
      Hit& hit,
      TFCSSimulationState& simulstate,
      const TFCSTruthState* truth,
      const TFCSExtrapolationState* extrapol) override;

  virtual bool operator==(const TFCSParametrizationBase& ref) const override;

  void Print(Option_t* option) const override;

protected:
  CaloGeo* m_geo;  //! do not persistify

private:
  ClassDefOverride(TFCSHitCellMapping, 1)  // TFCSHitCellMapping
};
}  // namespace FastCaloSim::Core

#endif
