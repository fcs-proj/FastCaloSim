/*
  Copyright (C) 2002-2020 CERN for the benefit of the ATLAS collaboration
*/

#ifndef TFCSHitCellMapping_h
#define TFCSHitCellMapping_h

#include "FastCaloSim/Core/TFCSLateralShapeParametrizationHitBase.h"

class ICaloGeometry;

class TFCSHitCellMapping : public TFCSLateralShapeParametrizationHitBase
{
public:
  TFCSHitCellMapping(const char* name = nullptr,
                     const char* title = nullptr,
                     ICaloGeometry* geo = nullptr);

  virtual void set_geometry(ICaloGeometry* geo) override { m_geo = geo; };
  ICaloGeometry* get_geometry() { return m_geo; };

  /// fills all hits into calorimeter cells
  virtual FCSReturnCode simulate_hit(
      Hit& hit,
      TFCSSimulationState& simulstate,
      const TFCSTruthState* truth,
      const TFCSExtrapolationState* extrapol) override;

  virtual bool operator==(const TFCSParametrizationBase& ref) const override;

  void Print(Option_t* option) const override;

protected:
  ICaloGeometry* m_geo;  //! do not persistify

private:
  ClassDefOverride(TFCSHitCellMapping, 1)  // TFCSHitCellMapping
};

#endif
