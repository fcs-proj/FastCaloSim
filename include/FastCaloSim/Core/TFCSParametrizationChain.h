// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#pragma once

#include "FastCaloSim/Core/TFCSParametrization.h"

class TFCSParametrizationChain : public TFCSParametrization
{
public:
  TFCSParametrizationChain(const char* name = nullptr,
                           const char* title = nullptr)
      : TFCSParametrization(name, title)
  {
  }
  TFCSParametrizationChain(const TFCSParametrizationChain& ref)
      : TFCSParametrization(ref.GetName(), ref.GetTitle())
      , m_chain(ref.chain())
  {
  }

  /// Status bit for chain persistency
  enum FCSSplitChainObjects
  {
    kSplitChainObjects =
        BIT(16),  ///< Set this bit in the TObject bit field if the
                  ///< TFCSParametrizationBase objects in the chain should be
                  ///< written as separate keys into the root file instead of
                  ///< directly writing the objects. This is needed if the sum
                  ///< of all objects in the chain use >1GB of memory, which
                  ///< can't be handled by TBuffer. Drawback is that identical
                  ///< objects will get stored as multiple instances
    kRetryChainFromStart = BIT(17)
  };

  auto SplitChainObjects() const -> bool { return TestBit(kSplitChainObjects); }
  void set_SplitChainObjects() { SetBit(kSplitChainObjects); }
  void reset_SplitChainObjects() { ResetBit(kSplitChainObjects); }

  auto RetryChainFromStart() const -> bool
  {
    return TestBit(kRetryChainFromStart);
  }
  void set_RetryChainFromStart() { SetBit(kRetryChainFromStart); }
  void reset_RetryChainFromStart() { ResetBit(kRetryChainFromStart); }

  using Chain_t = std::vector<TFCSParametrizationBase*>;
  auto size() const -> unsigned int override
  {
    return static_cast<unsigned int>(m_chain.size());
  }
  auto operator[](unsigned int ind) const
      -> const TFCSParametrizationBase* override
  {
    return m_chain[ind];
  }
  auto operator[](unsigned int ind) -> TFCSParametrizationBase* override
  {
    return m_chain[ind];
  }
  void set_daughter(unsigned int ind, TFCSParametrizationBase* param) override
  {
    m_chain.at(ind) = param;
  }
  auto chain() const -> const Chain_t& { return m_chain; }
  auto chain() -> Chain_t& { return m_chain; }
  void push_back(const Chain_t::value_type& param)
  {
    m_chain.push_back(param);
    recalc();
  }

  auto is_match_Ekin_bin(int Ekin_bin) const -> bool override;
  auto is_match_calosample(int calosample) const -> bool override;

  auto simulate(TFCSSimulationState& simulstate,
                const TFCSTruthState* truth,
                const TFCSExtrapolationState* extrapol) const
      -> FCSReturnCode override;

  void Print(Option_t* option = "") const override;

  // THIS CLASS HAS A CUSTOM STREAMER! CHANGES IN THE VERSIONING OR DATA TYPES
  // NEED TO BE IMPLEMENTED BY HAND! void
  // TFCSParametrizationChain::Streamer(TBuffer &R__b)
protected:
  void recalc_pdgid_intersect();
  void recalc_pdgid_union();

  void recalc_Ekin_intersect();
  void recalc_eta_intersect();
  void recalc_Ekin_eta_intersect();

  void recalc_Ekin_union();
  void recalc_eta_union();
  void recalc_Ekin_eta_union();

  /// Default is to call recalc_pdgid_intersect() and
  /// recalc_Ekin_eta_intersect()
  virtual void recalc();

  /// Simulate one parametrization of the chain. A parametrization asks to be
  /// simulated again by returning FCSRetry + n, which allows up to n retries.
  /// With RetryChainFromStart(), that request is returned to the caller
  /// instead, which then restarts the whole chain (simulate_with_restarts).
  auto simulate_and_retry(TFCSParametrizationBase* parametrization,
                          TFCSSimulationState& simulstate,
                          const TFCSTruthState* truth,
                          const TFCSExtrapolationState* extrapol) const
      -> FCSReturnCode;

  /// Run one pass over the chain, pass() -> FCSReturnCode, and repeat it from
  /// the start as long as it returns FCSRetry + n, at most n times. `name`
  /// identifies the caller in the log messages.
  template<typename Pass>
  auto simulate_with_restarts(const char* name, Pass pass) const
      -> FCSReturnCode;

  /// Retries from this one on are reported with a warning
  static auto first_warned_retry(int max_retries) -> int;

private:
  Chain_t m_chain;

  ClassDefOverride(TFCSParametrizationChain, 2)  // TFCSParametrizationChain
};

#include "FastCaloSim/Core/TFCSParametrizationChain.icc"
