/*
  Copyright (C) 2002-2020 CERN for the benefit of the ATLAS collaboration
*/

#ifndef TFCSHitCellMapping_h
#define TFCSHitCellMapping_h

#include "FastCaloSim/Core/TFCSLateralShapeParametrizationHitBase.h"
#include "FastCaloSim/Geometry/Cell.h"

class CaloGeo;

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

  auto mod_diff(double angle1, double angle2) -> double;
  auto dist(const Cell& cell, const Hit& hit) -> double;

  virtual bool operator==(const TFCSParametrizationBase& ref) const override;

  void Print(Option_t* option) const override;

protected:
  CaloGeo* m_geo;  //! do not persistify

private:
  ClassDefOverride(TFCSHitCellMapping, 1)  // TFCSHitCellMapping
};

#endif
