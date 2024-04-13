/*
  Copyright (C) 2002-2020 CERN for the benefit of the ATLAS collaboration
*/

#ifndef TFCSLateralShapeParametrizationHitChain_h
#define TFCSLateralShapeParametrizationHitChain_h

#include "FastCaloSim/TFCSLateralShapeParametrization.h"
#include "FastCaloSim/TFCSLateralShapeParametrizationHitBase.h"
#include <vector>

class ICaloGeometry;

class TFCSLateralShapeParametrizationHitChain
    : public TFCSLateralShapeParametrization {
public:
#ifdef USE_GPU
  // modify the construction function in order to get the geometry for FCS-GPU
  TFCSLateralShapeParametrizationHitChain(const char *name = nullptr,
                                          const char *title = nullptr,
                                          ICaloGeometry *geo = nullptr);
  // override the set_geometry function for FCS-GPU
  // the geometry is needed to assign the deposit energy to cell
  // will not compile by default
  virtual void set_geometry(ICaloGeometry *geo) override;
#else
  TFCSLateralShapeParametrizationHitChain(const char *name = nullptr,
                                          const char *title = nullptr);
#endif

  TFCSLateralShapeParametrizationHitChain(
      TFCSLateralShapeParametrizationHitBase *hitsim);

  virtual FCSReturnCode
  init_hit(TFCSLateralShapeParametrizationHitBase::Hit &hit,
           TFCSSimulationState &simulstate, const TFCSTruthState *truth,
           const TFCSExtrapolationState *extrapol) const;

  virtual FCSReturnCode
  simulate(TFCSSimulationState &simulstate, const TFCSTruthState *truth,
           const TFCSExtrapolationState *extrapol) const override;

  typedef std::vector<TFCSLateralShapeParametrizationHitBase *> Chain_t;
  virtual unsigned int size() const override;
  virtual const TFCSParametrizationBase *
  operator[](unsigned int ind) const override;
  virtual TFCSParametrizationBase *operator[](unsigned int ind) override;
  virtual void set_daughter(unsigned int ind,
                            TFCSParametrizationBase *param) override;
  const Chain_t &chain() const { return m_chain; };
  Chain_t &chain() { return m_chain; };
  void push_back(const Chain_t::value_type &value) {
    m_chain.push_back(value);
  };
  void push_back_init(const Chain_t::value_type &value);

  unsigned int get_nr_of_init() const { return m_ninit; };
  void set_nr_of_init(unsigned int ninit) { m_ninit = ninit; };

  /// set which instance should determine the number of hits
  virtual void
  set_number_of_hits_simul(TFCSLateralShapeParametrizationHitBase *sim) {
    m_number_of_hits_simul = sim;
  };

  /// Call get_number_of_hits() only once, as it could contain a random number
  virtual int get_number_of_hits(TFCSSimulationState &simulstate,
                                 const TFCSTruthState *truth,
                                 const TFCSExtrapolationState *extrapol) const;

  /// Get hit energy from layer energy and number of hits
  virtual float get_E_hit(TFCSSimulationState &simulstate,
                          const TFCSTruthState *truth,
                          const TFCSExtrapolationState *extrapol) const;

  /// Give the effective size sigma^2 of the fluctuations that should be
  /// generated.
  virtual float
  get_sigma2_fluctuation(TFCSSimulationState &simulstate,
                         const TFCSTruthState *truth,
                         const TFCSExtrapolationState *extrapol) const;

  /// Get minimum and maximum value of weight for hit energy reweighting
  virtual float getMinWeight() const;
  virtual float getMaxWeight() const;

  static constexpr float s_max_sigma2_fluctuation = 1000; //! Do not persistify!

  void Print(Option_t *option = "") const override;

protected:
  void PropagateMSGLevel(MSG::Level level) const;

  Chain_t m_chain;

private:
#ifdef USE_GPU
  // the FCS-GPU will use this geometry to get the correct cell
  // will not compile by default
  ICaloGeometry *m_geo;
#endif

  TFCSLateralShapeParametrizationHitBase *m_number_of_hits_simul;
  unsigned int m_ninit = 0;

  ClassDefOverride(TFCSLateralShapeParametrizationHitChain,
                   2) // TFCSLateralShapeParametrizationHitChain
};

#endif
