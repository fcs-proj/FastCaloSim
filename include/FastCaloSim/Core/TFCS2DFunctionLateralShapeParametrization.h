// Copyright (c) 2024 CERN for the benefit of the FastCaloSim project

#ifndef TFCS2DFunctionLateralShapeParametrization_h
#define TFCS2DFunctionLateralShapeParametrization_h

#include "FastCaloSim/Core/TFCS2DFunction.h"
#include "FastCaloSim/Core/TFCSLateralShapeParametrizationHitBase.h"
#include "FastCaloSim/Core/TFCSTruthState.h"

class TH2;

class TFCS2DFunctionLateralShapeParametrization
    : public TFCSLateralShapeParametrizationHitBase
{
public:
  TFCS2DFunctionLateralShapeParametrization(const char* name = nullptr,
                                            const char* title = nullptr);
  ~TFCS2DFunctionLateralShapeParametrization();

  /// Status bit for FCS needs
  enum FCSStatusBits
  {
    k_phi_symmetric =
        BIT(15)  ///< Set this bit to simulate phi symmetric histograms
  };

  bool is_phi_symmetric() const { return TestBit(k_phi_symmetric); };
  virtual void set_phi_symmetric() { SetBit(k_phi_symmetric); };
  virtual void reset_phi_symmetric() { ResetBit(k_phi_symmetric); };

  /// set the desired number of hits
  void set_number_of_hits(float nhits);

  float get_number_of_expected_hits() const { return m_nhits; };

  /// default for this class is to simulate get_number_of_expected_hits() hits,
  /// which gives fluctuations sigma^2=1/get_number_of_expected_hits()
  virtual double get_sigma2_fluctuation(
      TFCSSimulationState& simulstate,
      const TFCSTruthState* truth,
      const TFCSExtrapolationState* extrapol) const override;

  /// default for this class is to simulate get_number_of_expected_hits() hits
  virtual int get_number_of_hits(
      TFCSSimulationState& simulstate,
      const TFCSTruthState* truth,
      const TFCSExtrapolationState* extrapol) const override;

  /// simulated one hit position with weight that should be put into simulstate
  /// sometime later all hit weights should be resacled such that their final
  /// sum is simulstate->E(sample) someone also needs to map all hits into cells
  virtual FCSReturnCode simulate_hit(
      Hit& hit,
      TFCSSimulationState& simulstate,
      const TFCSTruthState* truth,
      const TFCSExtrapolationState* extrapol) override;

  /// Init from function
  bool Initialize(TFCS2DFunction* func, float nhits = -1);

  TFCS2DFunction* function() { return m_function; };
  const TFCS2DFunction* function() const { return m_function; };

  virtual void Print(Option_t* option = "") const override;

protected:
  /// Histogram to be used for the shape simulation
  TFCS2DFunction* m_function;
  float m_nhits;

private:
  ClassDefOverride(TFCS2DFunctionLateralShapeParametrization,
                   2)  // TFCS2DFunctionLateralShapeParametrization
};

#endif
