// Copyright (c) 2024 CERN for the benefit of the FastCaloSim project

#include <iostream>

#include "FastCaloSim/Core/TFCS1DFunctionHistogram.h"

#include "TFile.h"
#include "TMath.h"
using namespace std;

//=============================================
//======= TFCS1DFunctionHistogram =========
//=============================================

TFCS1DFunctionHistogram::TFCS1DFunctionHistogram(TH1* hist, double cut_maxdev)
{
  Initialize(hist, cut_maxdev);
}

void TFCS1DFunctionHistogram::Initialize(TH1* hist, double cut_maxdev)
{
  smart_rebin_loop(hist, cut_maxdev);
}

double* TFCS1DFunctionHistogram::histo_to_array(TH1* hist)
{
  TH1D* h_clone = (TH1D*)hist->Clone("h_clone");
  h_clone->Scale(1.0 / h_clone->Integral());

  double* histoVals = new double[h_clone->GetNbinsX()];
  histoVals[0] = h_clone->GetBinContent(1);
  for (int i = 1; i < h_clone->GetNbinsX(); i++) {
    histoVals[i] = histoVals[i - 1] + h_clone->GetBinContent(i + 1);
  }
  delete h_clone;
  return histoVals;
}

double TFCS1DFunctionHistogram::sample_from_histo(TH1* hist, double random)
{
  double* histoVals = histo_to_array(hist);
  double value = 0.0;
  int chosenBin =
      (int)TMath::BinarySearch(hist->GetNbinsX(), histoVals, random);
  value = hist->GetBinCenter(chosenBin + 2);

  // cleanup
  delete[] histoVals;

  return value;
}

double TFCS1DFunctionHistogram::sample_from_histovalues(double random)
{
  double value = 0.0;

  TH1* hist = vector_to_histo();
  hist->SetName("hist");
  double* histoVals = histo_to_array(hist);
  int chosenBin =
      (int)TMath::BinarySearch(hist->GetNbinsX(), histoVals, random);
  value = hist->GetBinCenter(chosenBin + 2);

  return value;
}

TH1* TFCS1DFunctionHistogram::vector_to_histo()
{
  double* bins = new double[m_HistoBorders.size()];
  for (unsigned int i = 0; i < m_HistoBorders.size(); i++)
    bins[i] = m_HistoBorders[i];

  TH1* h_out = new TH1D("h_out", "h_out", m_HistoBorders.size() - 1, bins);
  for (int b = 1; b <= h_out->GetNbinsX(); b++)
    h_out->SetBinContent(b, m_HistoContents[b - 1]);

  delete[] bins;

  return h_out;
}

void TFCS1DFunctionHistogram::smart_rebin_loop(TH1* hist, double cut_maxdev)
{
  m_HistoContents.clear();
  m_HistoBorders.clear();

  double change = get_change(hist)
      * 1.000001;  // increase slightly for comparison of floats

  double maxdev = -1;

  TH1D* h_input = (TH1D*)hist->Clone("h_input");
  TH1D* h_output = nullptr;

  int i = 0;
  while (1) {
    TH1D* h_out;
    if (i == 0) {
      h_out = (TH1D*)h_input->Clone("h_out");
    } else {
      h_out = smart_rebin(h_input);
      h_out->SetName("h_out");
    }

    maxdev = get_maxdev(hist, h_out);
    maxdev *= 100.0;

    if (i % 100 == 0)
      ATH_MSG_INFO("Iteration nr. " << i << " -----> change " << change
                                    << " bins " << h_out->GetNbinsX()
                                    << " -> maxdev=" << maxdev);

    if (maxdev < cut_maxdev && h_out->GetNbinsX() > 5 && i < 1000) {
      delete h_input;
      h_input = (TH1D*)h_out->Clone("h_input");
      change = get_change(h_input) * 1.000001;
      delete h_out;
      i++;
    } else {
      h_output = (TH1D*)h_input->Clone("h_output");
      delete h_out;
      break;
    }
  }

  ATH_MSG_INFO("Info: Rebinned histogram has " << h_output->GetNbinsX()
                                               << " bins.");

  // store:

  for (int b = 1; b <= h_output->GetNbinsX(); b++)
    m_HistoBorders.push_back((float)h_output->GetBinLowEdge(b));
  m_HistoBorders.push_back((float)h_output->GetXaxis()->GetXmax());

  for (int b = 1; b < h_output->GetNbinsX(); b++)
    m_HistoContents.push_back(h_output->GetBinContent(b));
  m_HistoContents.push_back(1);
}

double TFCS1DFunctionHistogram::get_maxdev(TH1* h_in, TH1D* h_out)
{
  double maxdev = 0;
  for (int i = 1; i <= h_in->GetNbinsX(); i++) {
    int bin = h_out->FindBin(h_in->GetBinCenter(i));
    double val = fabs(h_out->GetBinContent(bin) - h_in->GetBinContent(i));
    if (val > maxdev)
      maxdev = val;
  }
  return maxdev;
}

double TFCS1DFunctionHistogram::get_change(TH1* histo)
{
  // return the smallest change between 2 bin contents, but don't check the last
  // bin, because that one never gets merged
  double minchange = 100.0;
  for (int b = 2; b < histo->GetNbinsX(); b++) {
    double diff = histo->GetBinContent(b) - histo->GetBinContent(b - 1);
    if (diff < minchange && diff > 0)
      minchange = diff;
  }

  return minchange;
}

TH1D* TFCS1DFunctionHistogram::smart_rebin(TH1D* h_input)
{
  TH1D* h_out1 = (TH1D*)h_input->Clone("h_out1");

  // get the smallest change
  double change = get_change(h_out1) * 1.00001;

  vector<double> content;
  vector<double> binborder;

  binborder.push_back(h_out1->GetXaxis()->GetXmin());

  int merged = 0;
  int skip = 0;
  int secondlastbin_merge = 0;
  for (int b = 1; b < h_out1->GetNbinsX() - 1; b++)  // never touch the last bin
  {
    double thisBin = h_out1->GetBinContent(b);
    double nextBin = h_out1->GetBinContent(b + 1);
    double width = h_out1->GetBinWidth(b);
    double nextwidth = h_out1->GetBinWidth(b + 1);
    double diff = fabs(nextBin - thisBin);
    if (!skip && (diff > change || merged)) {
      binborder.push_back(h_out1->GetBinLowEdge(b + 1));
      content.push_back(thisBin);
    }
    skip = 0;
    if (diff <= change && !merged) {
      double sum = thisBin * width + nextBin * nextwidth;
      double sumwidth = width + nextwidth;
      binborder.push_back(h_out1->GetBinLowEdge(b + 2));
      content.push_back(sum / sumwidth);
      merged = 1;
      skip = 1;
      if (b == (h_out1->GetNbinsX() - 2))
        secondlastbin_merge = 1;
    }
  }  // for b
  if (!secondlastbin_merge) {
    binborder.push_back(h_out1->GetBinLowEdge(h_out1->GetNbinsX()));
    content.push_back(h_out1->GetBinContent(h_out1->GetNbinsX() - 1));
  }
  binborder.push_back(h_out1->GetXaxis()->GetXmax());
  content.push_back(h_out1->GetBinContent(h_out1->GetNbinsX()));

  double* bins = new double[content.size() + 1];
  for (unsigned int i = 0; i < binborder.size(); i++)
    bins[i] = binborder[i];

  TH1D* h_out2 = new TH1D("h_out2", "h_out2", content.size(), bins);
  for (unsigned int b = 1; b <= content.size(); b++)
    h_out2->SetBinContent(b, content[b - 1]);

  delete[] bins;
  delete h_out1;

  return h_out2;
}

double TFCS1DFunctionHistogram::rnd_to_fct(double rnd) const
{
  double value2 = get_inverse(rnd);

  return value2;
}

double TFCS1DFunctionHistogram::linear(
    double y1, double y2, double x1, double x2, double y)
{
  double x = -1;

  double eps = 0.0000000001;
  if ((y2 - y1) < eps)
    x = x1;
  else {
    double m = (y2 - y1) / (x2 - x1);
    double n = y1 - m * x1;
    x = (y - n) / m;
  }

  return x;
}

double TFCS1DFunctionHistogram::non_linear(
    double y1, double y2, double x1, double x2, double y)
{
  double x = -1;
  double eps = 0.0000000001;
  if ((y2 - y1) < eps)
    x = x1;
  else {
    double b = (x1 - x2) / (sqrt(y1) - sqrt(y2));
    double a = x1 - b * sqrt(y1);
    x = a + b * sqrt(y);
  }
  return x;
}

double TFCS1DFunctionHistogram::get_inverse(double rnd) const
{
  double value = 0.;

  if (rnd < m_HistoContents[0]) {
    double x1 = m_HistoBorders[0];
    double x2 = m_HistoBorders[1];
    double y1 = 0;
    double y2 = m_HistoContents[0];
    double x = non_linear(y1, y2, x1, x2, rnd);
    value = x;
  } else {
    // find the first HistoContent element that is larger than rnd:
    vector<float>::const_iterator larger_element =
        std::upper_bound(m_HistoContents.begin(), m_HistoContents.end(), rnd);
    size_t index = larger_element - m_HistoContents.begin();
    if (index >= m_HistoContents.size()) {
      --index;
    }
    double y = m_HistoContents[index];
    double x1 = m_HistoBorders[index];
    double x2 = m_HistoBorders[index + 1];
    double y1 = m_HistoContents[index - 1];
    double y2 = y;
    if ((index + 1) == (m_HistoContents.size() - 1)) {
      x2 = m_HistoBorders[m_HistoBorders.size() - 1];
      y2 = m_HistoContents[m_HistoContents.size() - 1];
    }
    double x = non_linear(y1, y2, x1, x2, rnd);
    value = x;
  }

  return value;
}
