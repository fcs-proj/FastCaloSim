// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#include <filesystem>
#include <memory>
#include <string>
#include <vector>

#include <TClass.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TParameter.h>
#include <gtest/gtest.h>
#include <unistd.h>

#include "FastCaloSim/Core/TFCSParametrizationBase.h"
#include "FastCaloSim/Core/TFCSPhiModulationCorrection.h"

namespace
{
// A correction created through ROOT, as when it is read from a
// parametrization file. The class is not exported, so the test owns it
// through the exported base class.
struct Correction
{
  std::unique_ptr<TFCSParametrizationBase> owner {
      static_cast<TFCSParametrizationBase*>(
          TClass::GetClass("TFCSPhiModulationCorrection")->New())};

  auto operator->() const -> TFCSPhiModulationCorrection*
  {
    return static_cast<TFCSPhiModulationCorrection*>(owner.get());
  }
};

// A file of its own for each test, as ctest runs the tests in parallel
auto file_path(const std::string& name) -> std::string
{
  return (std::filesystem::temp_directory_path()
          / (name + "_" + std::to_string(getpid()) + ".root"))
      .string();
}

const std::vector<float> weights {1.07F, 1.05F, 0.96F, 0.92F};

// The weights over phi within the cell, as in the phi modulation inputs
void write_th1(const std::string& path)
{
  TFile file(path.c_str(), "RECREATE");
  TH1F hist("hWt_Layer0", "", 4, 0.0, 0.1);
  hist.SetDirectory(nullptr);
  for (std::size_t i = 0; i < weights.size(); ++i) {
    hist.SetBinContent(static_cast<int>(i) + 1,
                       static_cast<double>(weights[i]));
  }
  hist.Write();
  TParameter<double> shift("energy_shift", 0.5);
  shift.Write();
}
}  // namespace

TEST(PhiModulationTest, LoadsOneDimensionalWeights)
{
  const std::string path = file_path("PhiModulationTest_TH1");
  write_th1(path);

  const Correction correction;
  correction->load_phi_modulation(path, 1, 0.2F, 0.25F);
  std::filesystem::remove(path);

  const auto& eta_mins = correction->get_min_eta();
  ASSERT_EQ(eta_mins.size(), 2U);
  EXPECT_EQ(eta_mins[1], (std::vector {0.2F, 0.25F}));

  const auto& modulation = correction->get_phi_modulation();
  ASSERT_EQ(modulation.size(), 2U);
  ASSERT_EQ(modulation[1].size(), 2U);
  EXPECT_EQ(modulation[1][0], weights);

  const auto& phi_mins = correction->get_min_phi();
  ASSERT_EQ(phi_mins[1][0].size(), weights.size() + 1);
  EXPECT_FLOAT_EQ(phi_mins[1][0].front(), 0.0F);
  EXPECT_FLOAT_EQ(phi_mins[1][0].back(), 0.1F);

  EXPECT_FLOAT_EQ(correction->get_energy_shift()[1][0], 0.5F);
}

// A file without a one-dimensional hWt_Layer0 is rejected and leaves the
// tables unchanged
TEST(PhiModulationTest, RejectsOtherHistograms)
{
  const std::string path = file_path("PhiModulationTest_TH2");
  {
    TFile file(path.c_str(), "RECREATE");
    TH2F hist("hWt_Layer0", "", 4, 0.0, 0.1, 1, 0.0, 1.0);
    hist.SetDirectory(nullptr);
    hist.Write();
  }

  const Correction correction;
  correction->load_phi_modulation(path, 1, 0.2F, 0.25F);
  correction->load_phi_modulation(
      file_path("PhiModulationTest_missing"), 1, 0.2F, 0.25F);
  std::filesystem::remove(path);

  EXPECT_TRUE(correction->get_min_eta().empty());
  EXPECT_TRUE(correction->get_phi_modulation().empty());
}
