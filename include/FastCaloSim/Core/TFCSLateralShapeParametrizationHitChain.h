// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#ifndef TFCSLateralShapeParametrizationHitChain_h
#define TFCSLateralShapeParametrizationHitChain_h

#include "FastCaloSim/Core/TFCSLateralShapeParametrization.h"
#include "FastCaloSim/Core/TFCSLateralShapeParametrizationHitBase.h"

class CaloGeo;

class TFCSLateralShapeParametrizationHitChain
    : public TFCSLateralShapeParametrization
{
public:
  TFCSLateralShapeParametrizationHitChain(const char* name = nullptr,
                                          const char* title = nullptr);

  TFCSLateralShapeParametrizationHitChain(
      TFCSLateralShapeParametrizationHitBase* hitsim);

  virtual auto init_hit(TFCSLateralShapeParametrizationHitBase::Hit& hit,
                        TFCSSimulationState& simulstate,
                        const TFCSTruthState* truth,
                        const TFCSExtrapolationState* extrapol) const
      -> FCSReturnCode;

  auto simulate(TFCSSimulationState& simulstate,
                const TFCSTruthState* truth,
                const TFCSExtrapolationState* extrapol) const
      -> FCSReturnCode override;

  using Chain_t = std::vector<TFCSLateralShapeParametrizationHitBase*>;
  auto size() const -> unsigned int override;
  auto operator[](unsigned int ind) const
      -> const TFCSParametrizationBase* override;
  auto operator[](unsigned int ind) -> TFCSParametrizationBase* override;
  void set_daughter(unsigned int ind, TFCSParametrizationBase* param) override;
  auto chain() const -> const Chain_t& { return m_chain; };
  auto chain() -> Chain_t& { return m_chain; };
  void push_back(const Chain_t::value_type& value)
  {
    m_chain.push_back(value);
  };
  void push_back_init(const Chain_t::value_type& value);

  auto get_nr_of_init() const -> unsigned int { return m_ninit; };
  void set_nr_of_init(unsigned int ninit) { m_ninit = ninit; };

  /// set which instance should determine the number of hits
  virtual void set_number_of_hits_simul(
      TFCSLateralShapeParametrizationHitBase* sim)
  {
    m_number_of_hits_simul = sim;
  };

  /// Call get_number_of_hits() only once, as it could contain a random number
  virtual auto get_number_of_hits(TFCSSimulationState& simulstate,
                                  const TFCSTruthState* truth,
                                  const TFCSExtrapolationState* extrapol) const
      -> int;

  /// Get hit energy from layer energy and number of hits
  virtual auto get_E_hit(TFCSSimulationState& simulstate,
                         const TFCSTruthState* truth,
                         const TFCSExtrapolationState* extrapol) const -> float;

  /// Give the effective size sigma^2 of the fluctuations that should be
  /// generated.
  virtual auto get_sigma2_fluctuation(
      TFCSSimulationState& simulstate,
      const TFCSTruthState* truth,
      const TFCSExtrapolationState* extrapol) const -> float;

  /// Get minimum and maximum value of weight for hit energy reweighting
  virtual auto getMinWeight() const -> float;
  virtual auto getMaxWeight() const -> float;

  static constexpr float s_max_sigma2_fluctuation =
      1000;  //! Do not persistify!

  void Print(Option_t* option = "") const override;

protected:
  void PropagateMSGLevel(FCS_MSG::Level level) const;

  virtual auto check_all_hits_simulated(
      TFCSLateralShapeParametrizationHitBase::Hit& hit,
      TFCSSimulationState& simulstate,
      const TFCSTruthState* truth,
      const TFCSExtrapolationState* extrapol,
      bool success) const -> bool;

  Chain_t m_chain;

private:
  TFCSLateralShapeParametrizationHitBase* m_number_of_hits_simul;
  unsigned int m_ninit = 0;

  ClassDefOverride(TFCSLateralShapeParametrizationHitChain,
                   2)  // TFCSLateralShapeParametrizationHitChain
};

#endif
