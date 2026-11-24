// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#include <memory>

#include <TClass.h>
#include <TMemFile.h>
#include <TObject.h>
#include <gtest/gtest.h>

#include "FastCaloSim/Core/TFCSParametrizationBase.h"
#include "FastCaloSim/Core/TFCSParametrizationChain.h"

namespace
{
// Created through ROOT, as when read from a parametrization file. The class is
// not exported (and its header needs libxml2), so the test owns it through the
// exported base class and asks its chain base class.
auto make_gan() -> std::unique_ptr<TFCSParametrizationBase>
{
  return std::unique_ptr<TFCSParametrizationBase>(
      static_cast<TFCSParametrizationBase*>(
          TClass::GetClass("TFCSEnergyAndHitGANV2")->New()));
}

auto restarts(const TFCSParametrizationBase& param) -> bool
{
  return static_cast<const TFCSParametrizationChain&>(param)
      .RetryChainFromStart();
}
}  // namespace

// A new GAN does not restart its chain
TEST(GANStreamerTest, NewGANDoesNotRestart)
{
  EXPECT_FALSE(restarts(*make_gan()));
}

// Older files have BIT(17), the bit of RetryChainFromStart(), set on every GAN
// for a flag that was never used: it is cleared when reading
TEST(GANStreamerTest, ReadingClearsTheOldFlag)
{
  const auto gan = make_gan();
  gan->SetBit(BIT(17));
  ASSERT_TRUE(restarts(*gan));

  TMemFile file("GANStreamerTest.root", "RECREATE");
  file.WriteObjectAny(gan.get(), "TFCSEnergyAndHitGANV2", "gan");
  const std::unique_ptr<TFCSParametrizationBase> read(
      file.Get<TFCSParametrizationBase>("gan"));
  ASSERT_NE(read, nullptr);
  EXPECT_FALSE(restarts(*read));
}
