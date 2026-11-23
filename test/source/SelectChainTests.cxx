// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#include <memory>

#include <TClass.h>
#include <gtest/gtest.h>

#include "FastCaloSim/Core/TFCSParametrizationBase.h"
#include "FastCaloSim/Core/TFCSParametrizationEkinSelectChain.h"

namespace
{
// A chain created through ROOT, as when it is read from a parametrization
// file. The chain classes are not exported, so the test owns the chain through
// the exported base class.
struct EkinChain
{
  std::unique_ptr<TFCSParametrizationBase> owner {
      static_cast<TFCSParametrizationBase*>(
          TClass::GetClass("TFCSParametrizationEkinSelectChain")->New())};

  auto operator->() const -> TFCSParametrizationEkinSelectChain*
  {
    return static_cast<TFCSParametrizationEkinSelectChain*>(owner.get());
  }
};
}  // namespace

// A second parametrization with the range of an existing bin goes into that
// bin, also for the last bin, whose upper edge is in no bin
TEST(SelectChainTest, PushBackInExistingBin)
{
  const EkinChain chain;
  TFCSParametrizationBase first("first", "first");
  TFCSParametrizationBase second("second", "second");
  TFCSParametrizationBase third("third", "third");
  TFCSParametrizationBase in_first_bin("in_first_bin", "in_first_bin");
  TFCSParametrizationBase in_last_bin("in_last_bin", "in_last_bin");

  EXPECT_EQ(chain->push_back_in_bin(&first, 1000.F, 2000.F), 0U);
  // Only one bin: it is the first and the last
  EXPECT_EQ(chain->push_back_in_bin(&second, 1000.F, 2000.F), 0U);
  EXPECT_EQ(chain->push_back_in_bin(&third, 2000.F, 4000.F), 1U);
  EXPECT_EQ(chain->push_back_in_bin(&in_first_bin, 1000.F, 2000.F), 0U);
  EXPECT_EQ(chain->push_back_in_bin(&in_last_bin, 2000.F, 4000.F), 1U);

  EXPECT_EQ(chain->get_number_of_bins(), 2U);
  EXPECT_EQ(chain->size(), 5U);
  EXPECT_EQ(chain->get_bin_low_edge(1), 2000.);
  EXPECT_EQ(chain->get_bin_up_edge(1), 4000.);

  // A range that splits the last bin is still rejected
  TFCSParametrizationBase splitting("splitting", "splitting");
  EXPECT_FALSE(chain->push_back_in_bin(&splitting, 2000.F, 3000.F));
  EXPECT_EQ(chain->size(), 5U);
}
