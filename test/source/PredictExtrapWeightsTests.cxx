// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#include <memory>

#include <TClass.h>
#include <gtest/gtest.h>

#include "FastCaloSim/Core/TFCSParametrizationBase.h"

namespace
{
// Created through ROOT, as when read from a parametrization file. The class is
// not exported, so the test owns it through the exported base class.
auto make_weights() -> std::unique_ptr<TFCSParametrizationBase>
{
  return std::unique_ptr<TFCSParametrizationBase>(
      static_cast<TFCSParametrizationBase*>(
          TClass::GetClass("TFCSPredictExtrapWeights")->New()));
}
}  // namespace

// Without network configurations (as after reading with freemem()), two
// otherwise equal objects are compared without dereferencing them and are not
// taken as duplicates
TEST(PredictExtrapWeightsTest, CompareWithoutNetworkConfiguration)
{
  const auto first = make_weights();
  const auto second = make_weights();
  EXPECT_FALSE(*first == *second);
  EXPECT_TRUE(*first == *first);

  const TFCSParametrizationBase other_class("other", "other");
  EXPECT_FALSE(*first == other_class);
}
