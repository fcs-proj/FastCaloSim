// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#include <algorithm>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <random>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "FastCaloSim/Definitions/GuideTable.h"

#include <TClass.h>
#include <TFile.h>
#include <TObject.h>
#include <gtest/gtest.h>

#include "FastCaloSim/Core/TFCS1DFunctionInt32Histogram.h"
#include "FastCaloSim/Core/TFCS2DFunctionHistogram.h"
#include "FastCaloSim/Core/TFCSHistoLateralShapeParametrization.h"
#include "FastCaloSim/Core/TFCSHitCellMappingWiggle.h"
#include "FastCaloSim/Core/TFCSParametrizationBase.h"

namespace
{
using FastCaloSim::GuideTable;
using Int32Content = TFCS1DFunctionInt32Histogram::HistoContent_t;
// max_value (not exported from the library)
constexpr Int32Content max_value = std::numeric_limits<Int32Content>::max();

// The searches that the guide table replaces
auto full_search(const std::vector<float>& cumulative, double rnd)
    -> std::size_t
{
  const auto it = std::ranges::upper_bound(cumulative,
                                           rnd,
                                           std::ranges::less {},
                                           [](float c) -> double
                                           { return static_cast<double>(c); });
  return static_cast<std::size_t>(it - cumulative.begin());
}

auto full_search(const std::vector<Int32Content>& cumulative, Int32Content rnd)
    -> std::size_t
{
  const auto it = std::ranges::upper_bound(cumulative, rnd);
  return static_cast<std::size_t>(it - cumulative.begin());
}

// Bin contents with a given fraction of empty bins
auto random_contents(std::mt19937_64& rng,
                     std::size_t nbins,
                     double empty_fraction) -> std::vector<double>
{
  std::uniform_real_distribution<double> flat(0, 1);
  std::vector<double> contents(nbins);
  for (auto& content : contents) {
    content = flat(rng) < empty_fraction ? 0 : flat(rng);
  }
  if (std::ranges::all_of(contents, [](double c) -> bool { return c <= 0; })) {
    contents[nbins / 2] = 1;
  }
  return contents;
}

// Normalized like TFCS2DFunctionHistogram::Initialize
auto cumulative_float(const std::vector<double>& contents) -> std::vector<float>
{
  std::vector<float> cumulative;
  cumulative.reserve(contents.size());
  double integral = 0;
  for (const double content : contents) {
    integral += content;
    cumulative.push_back(static_cast<float>(integral));
  }
  for (float& c : cumulative) {
    c = static_cast<float>(static_cast<double>(c) / integral);
  }
  return cumulative;
}

// Scaled like TFCS1DFunctionInt32Histogram::Initialize
auto cumulative_int32(const std::vector<double>& contents)
    -> std::vector<Int32Content>
{
  std::vector<double> sums;
  sums.reserve(contents.size());
  double integral = 0;
  for (const double content : contents) {
    integral += content;
    sums.push_back(integral);
  }
  std::vector<Int32Content> cumulative;
  cumulative.reserve(sums.size());
  for (const double sum : sums) {
    cumulative.push_back(
        static_cast<Int32Content>(max_value * (sum / integral)));
  }
  return cumulative;
}

// Uniform random numbers, the cumulative values and their neighbors, the
// bucket edges and their neighbors, the ends of the range and invalid values
auto probes(std::mt19937_64& rng,
            const std::vector<float>& cumulative,
            std::size_t nbuckets) -> std::vector<double>
{
  std::vector<double> rnds;
  rnds.reserve(10000 + (3 * (cumulative.size() + nbuckets + 1)) + 11);
  std::uniform_real_distribution<double> flat(0, 1);
  for (int i = 0; i < 10000; ++i) {
    rnds.push_back(flat(rng));
  }
  const auto add_with_neighbors = [&rnds](double u) -> void
  {
    rnds.push_back(u);
    rnds.push_back(std::nextafter(u, -1.0));
    rnds.push_back(std::nextafter(u, 2.0));
  };
  for (const float c : cumulative) {
    add_with_neighbors(static_cast<double>(c));
  }
  for (std::size_t k = 0; k <= nbuckets; ++k) {
    add_with_neighbors(static_cast<double>(k) / static_cast<double>(nbuckets));
  }
  for (const double u : {0.0,
                         -0.0,
                         std::numeric_limits<double>::denorm_min(),
                         -std::numeric_limits<double>::denorm_min(),
                         -1.0,
                         1.0,
                         1.5,
                         std::nextafter(1.0, 0.0),
                         std::numeric_limits<double>::infinity(),
                         -std::numeric_limits<double>::infinity(),
                         std::numeric_limits<double>::quiet_NaN()})
  {
    rnds.push_back(u);
  }
  return rnds;
}

auto probes(std::mt19937_64& rng,
            const std::vector<Int32Content>& cumulative,
            std::size_t nbuckets) -> std::vector<Int32Content>
{
  std::vector<Int32Content> rnds;
  rnds.reserve(10000 + (3 * (cumulative.size() + nbuckets + 2)));
  std::uniform_int_distribution<Int32Content> flat;
  for (int i = 0; i < 10000; ++i) {
    rnds.push_back(flat(rng));
  }
  const auto add_with_neighbors = [&rnds](std::uint64_t u) -> void
  {
    for (const std::uint64_t v : {u - 1, u, u + 1}) {
      if (v <= std::numeric_limits<Int32Content>::max()) {
        rnds.push_back(static_cast<Int32Content>(v));
      }
    }
  };
  for (const Int32Content c : cumulative) {
    add_with_neighbors(c);
  }
  for (std::size_t k = 0; k <= nbuckets; ++k) {
    add_with_neighbors((std::uint64_t {k} << 32U) / nbuckets);
  }
  add_with_neighbors(0);
  return rnds;
}

// Fixed seeds, so that the tests are reproducible
auto make_rng(std::uint64_t seed) -> std::mt19937_64
{
  return std::mt19937_64(seed);
}

const std::vector<std::size_t> bin_counts = {
    1, 2, 3, 5, 8, 17, 100, 1000, 4097};

template<typename T>
void expect_full_search_results(const std::vector<T>& cumulative,
                                std::mt19937_64& rng,
                                const std::string& label)
{
  GuideTable<T> guide;
  guide.build(cumulative);
  const std::size_t nbuckets = guide.empty() ? 1 : guide.size() - 1;
  std::size_t mismatches = 0;
  for (const auto u : probes(rng, cumulative, nbuckets)) {
    if (guide.upper_bound(cumulative, u) != full_search(cumulative, u)) {
      ++mismatches;
    }
  }
  EXPECT_EQ(mismatches, 0) << label;
}

// The histograms are not exported from the library: create them through ROOT,
// as when they are read, and own them through TObject
template<typename T>
auto make_function(const char* class_name) -> std::unique_ptr<TObject>
{
  return std::unique_ptr<TObject>(
      static_cast<T*>(TClass::GetClass(class_name)->New()));
}

// rnd_to_fct as before the guide table
void reference_rnd_to_fct(const TFCS2DFunctionHistogram& hist,
                          double& valuex,
                          double& valuey,
                          double rnd0,
                          double rnd1)
{
  const auto& contents = hist.get_HistoContents();
  const auto& borders_x = hist.get_HistoBordersx();
  const auto& borders_y = hist.get_HistoBordersy();
  const auto content = [&](std::size_t bin) -> double
  { return static_cast<double>(contents[bin]); };
  const auto border_x = [&](std::size_t bin) -> double
  { return static_cast<double>(borders_x[bin]); };
  const auto border_y = [&](std::size_t bin) -> double
  { return static_cast<double>(borders_y[bin]); };
  const std::size_t ibin =
      std::min(full_search(contents, rnd0), contents.size() - 1);
  const std::size_t nbinsx = borders_x.size() - 1;
  const std::size_t biny = ibin / nbinsx;
  const std::size_t binx = ibin - nbinsx * biny;
  const double basecont = ibin > 0 ? content(ibin - 1) : 0;
  const double dcont = content(ibin) - basecont;
  if (dcont > 0) {
    valuex = border_x(binx)
        + (border_x(binx + 1) - border_x(binx)) * (rnd0 - basecont) / dcont;
  } else {
    valuex = border_x(binx) + (border_x(binx + 1) - border_x(binx)) / 2;
  }
  valuey = border_y(biny) + (border_y(biny + 1) - border_y(biny)) * rnd1;
}

auto reference_rnd_to_fct(const TFCS1DFunctionInt32Histogram& hist, double rnd)
    -> double
{
  const auto& contents = hist.get_HistoContents();
  const auto& borders = hist.get_HistoBordersx();
  const auto int_rnd = static_cast<Int32Content>(max_value * rnd);
  const std::size_t ibin =
      std::min(full_search(contents, int_rnd), contents.size() - 1);
  const Int32Content basecont = ibin > 0 ? contents[ibin - 1] : 0;
  const auto low_edge = static_cast<double>(borders[ibin]);
  const auto bin_width = static_cast<double>(borders[ibin + 1]) - low_edge;
  const Int32Content dcont = contents[ibin] - basecont;
  if (dcont > 0) {
    return low_edge
        + bin_width * static_cast<double>(int_rnd - basecont)
        / static_cast<double>(dcont);
  }
  return low_edge + bin_width / 2;
}

// The same bin gives the same value; the tolerance only allows for a
// different contraction of the floating-point operations (FMA) in the library
// and in the test. A different bin gives a different value, except at the
// common edge of two bins.
auto same_value(double a, double b) -> bool
{
  return std::abs(a - b) <= 1e-12 * std::max(1.0, std::abs(b));
}

auto count_mismatches(const TFCS2DFunctionHistogram& hist, std::mt19937_64& rng)
    -> std::size_t
{
  const TFCS2DFunction& function = hist;
  std::size_t mismatches = 0;
  std::uniform_real_distribution<double> flat(0, 1);
  const std::size_t nbuckets =
      hist.guide_table().empty() ? 1 : hist.guide_table().size() - 1;
  for (const double rnd0 : probes(rng, hist.get_HistoContents(), nbuckets)) {
    if (std::isnan(rnd0) || rnd0 < 0 || rnd0 >= 1) {
      continue;
    }
    const double rnd1 = flat(rng);
    double x = 0;
    double y = 0;
    double ref_x = 0;
    double ref_y = 0;
    function.rnd_to_fct(x, y, rnd0, rnd1);
    reference_rnd_to_fct(hist, ref_x, ref_y, rnd0, rnd1);
    if (!same_value(x, ref_x) || !same_value(y, ref_y)) {
      ++mismatches;
    }
  }
  return mismatches;
}

auto count_mismatches(const TFCS1DFunctionInt32Histogram& hist,
                      std::mt19937_64& rng) -> std::size_t
{
  const TFCS1DFunction& function = hist;
  std::size_t mismatches = 0;
  std::uniform_real_distribution<double> flat(0, 1);
  std::vector<double> rnds;
  rnds.reserve(10000 + hist.get_HistoContents().size() + 2);
  for (int i = 0; i < 10000; ++i) {
    rnds.push_back(flat(rng));
  }
  // Random numbers that end up exactly on the cumulative values
  for (const Int32Content c : hist.get_HistoContents()) {
    rnds.push_back(static_cast<double>(c) / max_value);
  }
  rnds.push_back(0);
  rnds.push_back(std::nextafter(1.0, 0.0));
  for (const double rnd : rnds) {
    if (!same_value(function.rnd_to_fct(rnd), reference_rnd_to_fct(hist, rnd)))
    {
      ++mismatches;
    }
  }
  return mismatches;
}

// All objects of type T in a parametrization tree (shared objects once)
template<typename T>
void collect(const TFCSParametrizationBase* param,
             std::set<const TObject*>& visited,
             std::vector<const T*>& found,
             const char* class_name)
{
  if (param == nullptr || !visited.insert(param).second) {
    return;
  }
  if constexpr (std::same_as<T, TFCS2DFunctionHistogram>) {
    if (param->InheritsFrom("TFCSHistoLateralShapeParametrization")) {
      found.push_back(
          &static_cast<const TFCSHistoLateralShapeParametrization*>(param)
               ->histogram());
    }
  } else {
    if (param->InheritsFrom("TFCSHitCellMappingWiggle")) {
      for (const TFCS1DFunction* function :
           static_cast<const TFCSHitCellMappingWiggle*>(param)->get_functions())
      {
        if (function != nullptr && function->InheritsFrom(class_name)) {
          found.push_back(static_cast<const T*>(function));
        }
      }
    }
  }
  for (std::size_t i = 0; i < param->size(); ++i) {
    collect(param->operator[](i), visited, found, class_name);
  }
}
}  // namespace

TEST(GuideTableTests, FloatSameAsFullSearch)
{
  auto rng = make_rng(1);
  for (const std::size_t nbins : bin_counts) {
    for (const double empty_fraction : {0.0, 0.5, 0.95}) {
      expect_full_search_results(
          cumulative_float(random_contents(rng, nbins, empty_fraction)),
          rng,
          "nbins=" + std::to_string(nbins)
              + " empty=" + std::to_string(empty_fraction));
    }
  }
}

TEST(GuideTableTests, FloatEdgeCases)
{
  auto rng = make_rng(2);
  const float below_one = std::nextafter(1.0F, 0.0F);
  const float above_one = std::nextafter(1.0F, 2.0F);
  const std::vector<std::vector<float>> cases = {
      {1},
      {0, 0, 0, 1},
      {1, 1, 1, 1},
      {0.25F, 0.25F, 0.5F, 0.5F, 0.5F, 1, 1},
      {0, 0.5F, below_one},
      {0, 0.5F, above_one},
      {0.125F, 0.25F, 0.375F, 0.5F, 0.625F, 0.75F, 0.875F, 1},
      {1e-30F, 2e-30F, 1},
      {0.5F, 0.25F, 1},  // not sorted: no table, full search
      {0.5F, std::numeric_limits<float>::quiet_NaN(), 1},  // no table
      {std::numeric_limits<float>::quiet_NaN()},  // no table
  };
  for (const auto& cumulative : cases) {
    expect_full_search_results(
        cumulative, rng, "case of size " + std::to_string(cumulative.size()));
  }

  GuideTable<float> guide;
  guide.build(cases[3]);
  EXPECT_FALSE(guide.empty());
  guide.build(cases[8]);
  EXPECT_TRUE(guide.empty());
  guide.build(cases[9]);
  EXPECT_TRUE(guide.empty());
  guide.build(std::vector<float> {});
  EXPECT_TRUE(guide.empty());
  EXPECT_EQ(guide.upper_bound(std::vector<float> {}, 0.5), 0);

  // A cumulative distribution of another size than the table is searched in
  // full
  const std::vector<float> longer = {0.1F, 0.2F, 0.3F, 0.4F, 0.5F, 1};
  guide.build(cases[3]);
  for (const double u : {0.05, 0.25, 0.45, 0.99}) {
    EXPECT_EQ(guide.upper_bound(longer, u), full_search(longer, u));
  }
}

TEST(GuideTableTests, Int32SameAsFullSearch)
{
  auto rng = make_rng(3);
  for (const std::size_t nbins : bin_counts) {
    for (const double empty_fraction : {0.0, 0.5, 0.95}) {
      expect_full_search_results(
          cumulative_int32(random_contents(rng, nbins, empty_fraction)),
          rng,
          "nbins=" + std::to_string(nbins)
              + " empty=" + std::to_string(empty_fraction));
    }
  }
}

TEST(GuideTableTests, Int32EdgeCases)
{
  auto rng = make_rng(4);
  constexpr Int32Content max = std::numeric_limits<Int32Content>::max();
  const std::vector<std::vector<Int32Content>> cases = {
      {max},
      {0, 0, 0, max},
      {max, max, max},
      {1, 1, max / 2, max / 2, max - 1},
      {0x40000000U, 0x80000000U, 0xC0000000U, max},
      {5, 3, max},  // not sorted: no table, full search
  };
  for (const auto& cumulative : cases) {
    expect_full_search_results(
        cumulative, rng, "case of size " + std::to_string(cumulative.size()));
  }
  GuideTable<Int32Content> guide;
  guide.build(cases[5]);
  EXPECT_TRUE(guide.empty());
}

TEST(GuideTableTests, Histogram2DSameAsBefore)
{
  auto rng = make_rng(5);
  // (bins in x, bins in y)
  const std::vector<std::pair<std::size_t, std::size_t>> shapes = {
      {1, 1}, {7, 1}, {20, 20}};
  for (const auto& [nbinsx, nbinsy] : shapes) {
    const std::size_t nbins = nbinsx * nbinsy;
    const auto owner =
        make_function<TFCS2DFunctionHistogram>("TFCS2DFunctionHistogram");
    auto& hist = *static_cast<TFCS2DFunctionHistogram*>(owner.get());
    for (std::size_t i = 0; i <= nbinsx; ++i) {
      hist.get_HistoBordersx().push_back(static_cast<float>(i) * 0.5F);
    }
    for (std::size_t i = 0; i <= nbinsy; ++i) {
      hist.get_HistoBordersy().push_back(static_cast<float>(i) * 0.25F);
    }
    hist.set_HistoContents(cumulative_float(random_contents(rng, nbins, 0.3)));
    EXPECT_FALSE(hist.guide_table().empty());
    EXPECT_EQ(count_mismatches(hist, rng), 0) << nbins << " bins";

    // A copy read back (TObject::Clone streams it) has its table again
    const std::unique_ptr<TObject> copy(hist.Clone());
    const auto& read = *static_cast<const TFCS2DFunctionHistogram*>(copy.get());
    EXPECT_EQ(read.get_HistoContents(), hist.get_HistoContents());
    EXPECT_EQ(read.guide_table().size(), hist.guide_table().size());
    EXPECT_EQ(count_mismatches(read, rng), 0) << nbins << " bins";
  }
}

TEST(GuideTableTests, Histogram1DInt32SameAsBefore)
{
  auto rng = make_rng(6);
  for (const std::size_t nbins : std::vector<std::size_t> {1, 9, 300}) {
    const auto owner = make_function<TFCS1DFunctionInt32Histogram>(
        "TFCS1DFunctionInt32Histogram");
    auto& hist = *static_cast<TFCS1DFunctionInt32Histogram*>(owner.get());
    for (std::size_t i = 0; i <= nbins; ++i) {
      hist.get_HistoBordersx().push_back(static_cast<float>(i) * 0.1F - 1);
    }
    hist.set_HistoContents(cumulative_int32(random_contents(rng, nbins, 0.3)));
    EXPECT_FALSE(hist.guide_table().empty());
    EXPECT_EQ(count_mismatches(hist, rng), 0) << nbins << " bins";

    const std::unique_ptr<TObject> copy(hist.Clone());
    const auto& read =
        *static_cast<const TFCS1DFunctionInt32Histogram*>(copy.get());
    EXPECT_EQ(read.get_HistoContents(), hist.get_HistoContents());
    EXPECT_EQ(read.guide_table().size(), hist.guide_table().size());
    EXPECT_EQ(count_mismatches(read, rng), 0) << nbins << " bins";
  }
}

// The NaN fix of set_geometry replaces the contents and rebuilds the table
TEST(GuideTableTests, LateralShapeNanFixRebuildsTable)
{
  const std::unique_ptr<TFCSParametrizationBase> owner(
      static_cast<TFCSHistoLateralShapeParametrization*>(
          TClass::GetClass("TFCSHistoLateralShapeParametrization")->New()));
  auto& shape =
      *static_cast<TFCSHistoLateralShapeParametrization*>(owner.get());
  auto& hist = shape.histogram();
  hist.get_HistoBordersx() = {0, 1, 2, 3};
  hist.get_HistoBordersy() = {0, 1};
  const float nan = std::numeric_limits<float>::quiet_NaN();
  hist.set_HistoContents({0.25F, nan, nan});
  EXPECT_TRUE(hist.guide_table().empty());

  owner->set_geometry(nullptr);
  EXPECT_EQ(hist.get_HistoContents(), (std::vector<float> {0.25F, 1, 1}));
  EXPECT_FALSE(hist.guide_table().empty());
  auto rng = make_rng(7);
  EXPECT_EQ(count_mismatches(hist, rng), 0);
}

// The NaN fix keeps the cumulative distribution non-decreasing: probability 0
// before the first valid bin, 1 after the last, and none added in between
TEST(GuideTableTests, LateralShapeNanFixKeepsDistributionNonDecreasing)
{
  const float nan = std::numeric_limits<float>::quiet_NaN();
  const std::vector<std::pair<std::vector<float>, std::vector<float>>> cases {
      {{nan, 0.2F, 0.8F, nan}, {0, 0.2F, 0.8F, 1}},
      {{nan, nan, 0.5F, 1}, {0, 0, 0.5F, 1}},
      {{0.2F, nan, 0.6F, nan}, {0.2F, 0.2F, 0.6F, 1}},
      {{nan, nan, nan, nan}, {1, 1, 1, 1}},
  };
  for (const auto& [input, expected] : cases) {
    const std::unique_ptr<TFCSParametrizationBase> owner(
        static_cast<TFCSHistoLateralShapeParametrization*>(
            TClass::GetClass("TFCSHistoLateralShapeParametrization")->New()));
    auto& hist = static_cast<TFCSHistoLateralShapeParametrization*>(owner.get())
                     ->histogram();
    hist.get_HistoBordersx() = {0, 1, 2, 3, 4};
    hist.get_HistoBordersy() = {0, 1};
    hist.set_HistoContents(input);

    owner->set_geometry(nullptr);
    EXPECT_EQ(hist.get_HistoContents(), expected);
    EXPECT_TRUE(std::ranges::is_sorted(hist.get_HistoContents()));
  }
}

// Every histogram read from a parametrization file has its table
TEST(GuideTableTests, TablesBuiltWhenReadingParamFile)
{
  const std::unique_ptr<TFile> file(
      TFile::Open(TEST_BASE_DIR
                  "data/param/barrel/"
                  "TFCSParam_pid22_pid211_Mom16384_131072_eta_15_30.root"));
  ASSERT_NE(file, nullptr);
  const auto* param = file->Get<TFCSParametrizationBase>("SelPDGID");
  ASSERT_NE(param, nullptr);

  auto rng = make_rng(8);
  std::set<const TObject*> visited;
  std::vector<const TFCS2DFunctionHistogram*> hists_2d;
  collect(param, visited, hists_2d, "TFCS2DFunctionHistogram");
  EXPECT_FALSE(hists_2d.empty());
  for (const auto* hist : hists_2d) {
    EXPECT_FALSE(hist->get_HistoContents().empty());
    EXPECT_FALSE(hist->guide_table().empty());
    EXPECT_EQ(count_mismatches(*hist, rng), 0);
  }

  visited.clear();
  std::vector<const TFCS1DFunctionInt32Histogram*> hists_1d;
  collect(param, visited, hists_1d, "TFCS1DFunctionInt32Histogram");
  EXPECT_FALSE(hists_1d.empty());
  for (const auto* hist : hists_1d) {
    EXPECT_FALSE(hist->guide_table().empty());
    EXPECT_EQ(count_mismatches(*hist, rng), 0);
  }
}
