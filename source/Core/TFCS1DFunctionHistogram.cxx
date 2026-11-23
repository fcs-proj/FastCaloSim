// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#include <algorithm>
#include <cmath>
#include <vector>

#include "FastCaloSim/Core/TFCS1DFunctionHistogram.h"

#include "FastCaloSim/Definitions/Index.h"

TFCS1DFunctionHistogram::TFCS1DFunctionHistogram(TH1* hist, double cut_maxdev)
{
  Initialize(hist, cut_maxdev);
}

void TFCS1DFunctionHistogram::Initialize(TH1* hist, double cut_maxdev)
{
  smart_rebin_loop(hist, cut_maxdev);
}

void TFCS1DFunctionHistogram::smart_rebin_loop(TH1* hist, double cut_maxdev)
{
  m_HistoContents.clear();
  m_HistoBorders.clear();

  double change = get_change(hist)
      * 1.000001;  // increase slightly for comparison of floats

  double maxdev = -1;

  std::unique_ptr<TH1D> h_input(static_cast<TH1D*>(hist->Clone("h_input")));

  int i = 0;
  while (true) {
    std::unique_ptr<TH1D> h_out;
    if (i == 0) {
      h_out.reset(static_cast<TH1D*>(h_input->Clone("h_out")));
    } else {
      h_out = smart_rebin(h_input.get());
      h_out->SetName("h_out");
    }

    maxdev = get_maxdev(hist, h_out.get());
    maxdev *= 100.0;

    if (i % 100 == 0) {
      FCS_MSG_INFO("Iteration nr. " << i << " -----> change " << change
                                    << " bins " << h_out->GetNbinsX()
                                    << " -> maxdev=" << maxdev);
    }

    if (maxdev < cut_maxdev && h_out->GetNbinsX() > 5 && i < 1000) {
      h_input = std::move(h_out);
      h_input->SetName("h_input");
      change = get_change(h_input.get()) * 1.000001;
      i++;
    } else {
      break;
    }
  }

  FCS_MSG_INFO("Info: Rebinned histogram has " << h_input->GetNbinsX()
                                               << " bins.");

  // store:

  for (int b = 1; b <= h_input->GetNbinsX(); b++) {
    m_HistoBorders.push_back(static_cast<float>(h_input->GetBinLowEdge(b)));
  }
  m_HistoBorders.push_back(static_cast<float>(h_input->GetXaxis()->GetXmax()));

  for (int b = 1; b < h_input->GetNbinsX(); b++) {
    m_HistoContents.push_back(static_cast<float>(h_input->GetBinContent(b)));
  }
  m_HistoContents.push_back(1);
}

auto TFCS1DFunctionHistogram::get_maxdev(TH1* h_in, TH1D* h_out) -> double
{
  double maxdev = 0;
  for (int i = 1; i <= h_in->GetNbinsX(); i++) {
    int const bin = h_out->FindBin(h_in->GetBinCenter(i));
    double const val =
        std::fabs(h_out->GetBinContent(bin) - h_in->GetBinContent(i));
    maxdev = std::max(val, maxdev);
  }
  return maxdev;
}

auto TFCS1DFunctionHistogram::get_change(TH1* histo) -> double
{
  // return the smallest change between 2 bin contents, but don't check the last
  // bin, because that one never gets merged
  double minchange = 100.0;
  for (int b = 2; b < histo->GetNbinsX(); b++) {
    double const diff = histo->GetBinContent(b) - histo->GetBinContent(b - 1);
    if (diff < minchange && diff > 0) {
      minchange = diff;
    }
  }

  return minchange;
}

auto TFCS1DFunctionHistogram::smart_rebin(TH1D* h_input)
    -> std::unique_ptr<TH1D>
{
  std::unique_ptr<TH1D> h_out1(static_cast<TH1D*>(h_input->Clone("h_out1")));

  // get the smallest change
  double const change = get_change(h_out1.get()) * 1.00001;

  std::vector<double> content;
  std::vector<double> binborder;

  binborder.push_back(h_out1->GetXaxis()->GetXmin());

  bool merged = false;
  bool skip = false;
  bool secondlastbin_merge = false;
  for (int b = 1; b < h_out1->GetNbinsX() - 1; b++)  // never touch the last bin
  {
    double const thisBin = h_out1->GetBinContent(b);
    double const nextBin = h_out1->GetBinContent(b + 1);
    double const width = h_out1->GetBinWidth(b);
    double const nextwidth = h_out1->GetBinWidth(b + 1);
    double const diff = std::fabs(nextBin - thisBin);
    if (!skip && (diff > change || merged)) {
      binborder.push_back(h_out1->GetBinLowEdge(b + 1));
      content.push_back(thisBin);
    }
    skip = false;
    if (diff <= change && !merged) {
      double const sum = thisBin * width + nextBin * nextwidth;
      double const sumwidth = width + nextwidth;
      binborder.push_back(h_out1->GetBinLowEdge(b + 2));
      content.push_back(sum / sumwidth);
      merged = true;
      skip = true;
      if (b == (h_out1->GetNbinsX() - 2)) {
        secondlastbin_merge = true;
      }
    }
  }  // for b
  if (!secondlastbin_merge) {
    binborder.push_back(h_out1->GetBinLowEdge(h_out1->GetNbinsX()));
    content.push_back(h_out1->GetBinContent(h_out1->GetNbinsX() - 1));
  }
  binborder.push_back(h_out1->GetXaxis()->GetXmax());
  content.push_back(h_out1->GetBinContent(h_out1->GetNbinsX()));

  auto h_out2 = std::make_unique<TH1D>(
      "h_out2", "h_out2", content.size(), binborder.data());
  for (unsigned int b = 1; b <= content.size(); b++) {
    h_out2->SetBinContent(static_cast<Int_t>(b), content[b - 1]);
  }

  return h_out2;
}

auto TFCS1DFunctionHistogram::rnd_to_fct(double rnd) const -> double
{
  return get_inverse(rnd);
}

auto TFCS1DFunctionHistogram::non_linear(
    double y1, double y2, double x1, double x2, double y) -> double
{
  double x = -1;
  double const eps = 0.0000000001;
  if ((y2 - y1) < eps) {
    x = x1;
  } else {
    double const b = (x1 - x2) / (std::sqrt(y1) - std::sqrt(y2));
    double const a = x1 - b * std::sqrt(y1);
    x = a + b * std::sqrt(y);
  }
  return x;
}

auto TFCS1DFunctionHistogram::get_inverse(double rnd) const -> double
{
  // The histogram is stored in single precision
  const auto content = [this](std::size_t i) -> double
  { return static_cast<double>(m_HistoContents[i]); };
  const auto border = [this](std::size_t i) -> double
  { return static_cast<double>(m_HistoBorders[i]); };

  double value = 0.;

  if (rnd < content(0)) {
    double const x1 = border(0);
    double const x2 = border(1);
    double const y1 = 0;
    double const y2 = content(0);
    double const x = non_linear(y1, y2, x1, x2, rnd);
    value = x;
  } else {
    // find the first HistoContent element that is larger than rnd:
    auto larger_element = std::ranges::upper_bound(m_HistoContents, rnd);
    auto index = FastCaloSim::index_of(m_HistoContents, larger_element);
    if (index >= m_HistoContents.size()) {
      --index;
    }
    // Interpolate within the bin of rnd, between its edges and the
    // cumulative probabilities at them
    double const x1 = border(index);
    double const x2 = border(index + 1);
    double const y1 = content(index - 1);
    double const y2 = content(index);
    double const x = non_linear(y1, y2, x1, x2, rnd);
    value = x;
  }

  return value;
}
