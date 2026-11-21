// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#include <memory>
#include <vector>

#include <TClass.h>
#include <gtest/gtest.h>

#include "FastCaloSim/Core/TFCSParametrizationBase.h"
#include "FastCaloSim/Core/TFCSParametrizationBinnedChain.h"
#include "FastCaloSim/Core/TFCSParametrizationChain.h"
#include "FastCaloSim/Core/TFCSSimulationState.h"
#include "FastCaloSim/Core/TFCSTruthState.h"

namespace
{
// Returns the scripted status codes one after the other, then FCSSuccess
class ScriptedParametrization : public TFCSParametrizationBase
{
public:
  explicit ScriptedParametrization(std::vector<int> script)
      : m_script(std::move(script))
  {
  }

  auto simulate(TFCSSimulationState& /*simulstate*/,
                const TFCSTruthState* /*truth*/,
                const TFCSExtrapolationState* /*extrapol*/) const
      -> FCSReturnCode override
  {
    const int code = m_calls < m_script.size() ? m_script[m_calls] : FCSSuccess;
    ++m_calls;
    return static_cast<FCSReturnCode>(code);
  }

  auto calls() const -> std::size_t { return m_calls; }

private:
  std::vector<int> m_script;
  mutable std::size_t m_calls = 0;
};

// The chain classes are created through ROOT, as when they are read from a
// parametrization file
template<typename T>
auto make_chain(const char* class_name) -> std::unique_ptr<T>
{
  return std::unique_ptr<T>(
      static_cast<T*>(TClass::GetClass(class_name)->New()));
}

auto simulate(TFCSParametrizationBase& chain) -> FCSReturnCode
{
  TFCSSimulationState simulstate;
  simulstate.set_Ebin(1);
  TFCSTruthState truth;
  truth.set_pdgid(211);
  return chain.simulate(simulstate, &truth, nullptr);
}

constexpr auto retry(int n) -> int
{
  return FCSRetry + n;
}
}  // namespace

// A parametrization asking for FCSRetry + n is simulated again, alone
TEST(ChainRetryTest, ParametrizationRetried)
{
  auto chain = make_chain<TFCSParametrizationChain>("TFCSParametrizationChain");
  ScriptedParametrization first({});
  ScriptedParametrization second({retry(2), retry(2)});
  chain->push_back(&first);
  chain->push_back(&second);
  EXPECT_EQ(simulate(*chain), FCSSuccess);
  EXPECT_EQ(first.calls(), 1U);
  EXPECT_EQ(second.calls(), 3U);
}

TEST(ChainRetryTest, RetriesExhausted)
{
  auto chain = make_chain<TFCSParametrizationChain>("TFCSParametrizationChain");
  ScriptedParametrization param({retry(1), retry(1), retry(1)});
  chain->push_back(&param);
  EXPECT_EQ(simulate(*chain), FCSFatal);
  EXPECT_EQ(param.calls(), 2U);
}

TEST(ChainRetryTest, FatalStopsTheChain)
{
  auto chain = make_chain<TFCSParametrizationChain>("TFCSParametrizationChain");
  ScriptedParametrization first({FCSFatal});
  ScriptedParametrization second({});
  chain->push_back(&first);
  chain->push_back(&second);
  EXPECT_EQ(simulate(*chain), FCSFatal);
  EXPECT_EQ(second.calls(), 0U);
}

// With RetryChainFromStart, the whole chain is simulated again
TEST(ChainRetryTest, RestartFromStart)
{
  auto chain = make_chain<TFCSParametrizationChain>("TFCSParametrizationChain");
  chain->set_RetryChainFromStart();
  ScriptedParametrization first({});
  ScriptedParametrization second({retry(2)});
  ScriptedParametrization third({});
  chain->push_back(&first);
  chain->push_back(&second);
  chain->push_back(&third);
  EXPECT_EQ(simulate(*chain), FCSSuccess);
  EXPECT_EQ(first.calls(), 2U);
  EXPECT_EQ(second.calls(), 2U);
  EXPECT_EQ(third.calls(), 1U);
}

// Binned chain without parametrizations before the bins: the restart must run
// the parametrization of the bin again
TEST(ChainRetryTest, BinnedRestartWithoutParametrizationsBeforeBins)
{
  auto chain = make_chain<TFCSParametrizationBinnedChain>(
      "TFCSParametrizationEbinChain");
  chain->set_RetryChainFromStart();
  ScriptedParametrization in_bin({retry(2)});
  chain->push_back_in_bin(&in_bin, 1);
  EXPECT_EQ(simulate(*chain), FCSSuccess);
  EXPECT_EQ(in_bin.calls(), 2U);
}
