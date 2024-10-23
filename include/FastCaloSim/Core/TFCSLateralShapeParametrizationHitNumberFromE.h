// Copyright (c) 2024 CERN for the benefit of the FastCaloSim project

#ifndef TFCSLateralShapeParametrizationHitNumberFromE_h
#define TFCSLateralShapeParametrizationHitNumberFromE_h

#include "FastCaloSim/Core/TFCSLateralShapeParametrizationHitBase.h"
#include "TH2.h"

class TFCSLateralShapeParametrizationHitNumberFromE
    : public TFCSLateralShapeParametrizationHitBase
{
public:
  /// LAr: 10.1%/sqrt(E)
  ///    stochastic=0.101;
  ///    constant=0.002;
  /// HadEC: 21.4%/sqrt(E)
  ///    stochastic=0.214;
  ///    constant=0.0;
  /// TileCal: 56.4%/sqrt(E)
  ///    stochastic=0.564;
  ///    constant=0.055;
  /// FCAL:    28.5%/sqrt(E)
  ///    stochastic=0.285;
  ///    constant=0.035;
  TFCSLateralShapeParametrizationHitNumberFromE(const char* name = nullptr,
                                                const char* title = nullptr,
                                                double stochastic = 0.1,
                                                double constant = 0);

  TFCSLateralShapeParametrizationHitNumberFromE(const char* name,
                                                const char* title,
                                                double stochastic,
                                                double stochastic_hadron,
                                                double constant);

  /// Give the effective size sigma^2 of the fluctuations from the stochastic
  /// and constant term
  double get_sigma2_fluctuation(
      TFCSSimulationState& simulstate,
      const TFCSTruthState* truth,
      const TFCSExtrapolationState* extrapol) const override;

  int get_number_of_hits(TFCSSimulationState& simulstate,
                         const TFCSTruthState* truth,
                         const TFCSExtrapolationState* extrapol) const override;

  virtual bool operator==(const TFCSParametrizationBase& ref) const override;

  void Print(Option_t* option = "") const override;

protected:
  bool compare(const TFCSParametrizationBase& ref) const;

private:
  // Information for the fluctuation terms
  // The variation is calculated as:
  // sigma^2=[m_stochastic/sqrt(E/GeV)]^2 + [m_constant +
  // m_stochastic_hadron/sqrt(E/GeV)]^2
  double m_stochastic;
  double m_stochastic_hadron;
  double m_constant;

  ClassDefOverride(TFCSLateralShapeParametrizationHitNumberFromE,
                   2)  // TFCSLateralShapeParametrizationHitNumberFromE
};

#endif
