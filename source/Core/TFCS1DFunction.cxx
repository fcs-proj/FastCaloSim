// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-equal"

#include <cmath>
#include <cstddef>
#include <iostream>
#include <memory>

#include "FastCaloSim/Core/TFCS1DFunction.h"

#include "TH1.h"
#include "TRandom.h"

//=============================================
//======= TFCS1DFunction =========
//=============================================

void TFCS1DFunction::rnd_to_fct(float value[], const float rnd[]) const
{
  value[0] = static_cast<float>(rnd_to_fct(static_cast<double>(rnd[0])));
}

auto TFCS1DFunction::get_maxdev(TH1* h_input1, TH1* h_approx1) -> double
{
  std::unique_ptr<TH1D> h_input(static_cast<TH1D*>(h_input1->Clone("h_input")));
  std::unique_ptr<TH1D> h_approx(
      static_cast<TH1D*>(h_approx1->Clone("h_approx")));

  double maxdev = 0.0;

  // normalize the histos to the same area:
  double integral_input = h_input->Integral();
  double integral_approx = 0.0;
  for (int b = 1; b <= h_input->GetNbinsX(); b++) {
    integral_approx +=
        h_approx->GetBinContent(h_approx->FindBin(h_input->GetBinCenter(b)));
  }
  if (integral_approx == 0.0) {
    return 0.;
  }
  h_approx->Scale(integral_input / integral_approx);

  double ymax = h_approx->GetBinContent(h_approx->GetNbinsX())
      - h_approx->GetBinContent(h_approx->GetMinimumBin());
  for (int i = 1; i <= h_input->GetNbinsX(); i++) {
    double val = std::abs(h_approx->GetBinContent(
                              h_approx->FindBin(h_input->GetBinCenter(i)))
                          - h_input->GetBinContent(i))
        / ymax;
    if (val > maxdev)
      maxdev = val;
  }

  return maxdev * 100.0;
}

auto TFCS1DFunction::CheckAndIntegrate1DHistogram(
    const TH1* hist, std::vector<double>& integral_vec, int& first, int& last)
    -> double
{
  FastCaloSim::MLogging logger;
  Int_t nbins = hist->GetNbinsX();

  double integral = 0;
  integral_vec.resize(static_cast<std::size_t>(nbins));
  for (int ix = 1; ix <= nbins; ix++) {
    double binval = hist->GetBinContent(ix);
    if (binval < 0) {
      // Can't work if a bin is negative, forcing bins to 0 in this case
      double fraction = binval / hist->Integral();
      if (std::abs(fraction) > 1e-5) {
        FCS_MSG_NOCLASS(logger,
                        "Warning: bin content is negative in histogram "
                            << hist->GetName() << " : " << hist->GetTitle()
                            << " binval=" << binval << " " << fraction * 100
                            << "% of integral=" << hist->Integral()
                            << ". Forcing bin to 0.");
      }
      binval = 0;
    }
    integral += binval;
    integral_vec[static_cast<std::size_t>(ix - 1)] = integral;
  }

  for (first = 0; first < nbins; first++)
    if (integral_vec[static_cast<std::size_t>(first)] != 0)
      break;
  for (last = nbins - 1; last > 0; last--)
    if (integral_vec[static_cast<std::size_t>(last)] != integral)
      break;
  last++;

  if (integral <= 0) {
    FCS_MSG_NOCLASS(logger,
                    "Error: histogram "
                        << hist->GetName() << " : " << hist->GetTitle()
                        << " integral=" << integral << " is <=0");
  }
  return integral;
}

auto TFCS1DFunction::generate_histogram_random_slope(
    int nbinsx, double xmin, double xmax, double zerothreshold) -> TH1*
{
  TH1* hist = new TH1D("test_slope1D", "test_slope1D", nbinsx, xmin, xmax);
  hist->Sumw2();
  for (int ix = 1; ix <= nbinsx; ++ix) {
    double val = (0.5 + gRandom->Rndm()) * (nbinsx + ix);
    if (gRandom->Rndm() < zerothreshold)
      val = 0;
    hist->SetBinContent(ix, val);
    hist->SetBinError(ix, 0);
  }
  return hist;
}

auto TFCS1DFunction::generate_histogram_random_gauss(
    int nbinsx, int ntoy, double xmin, double xmax, double xpeak, double sigma)
    -> TH1*
{
  TH1* hist = new TH1D("test_gauss1D", "test_gauss1D", nbinsx, xmin, xmax);
  hist->Sumw2();
  for (int i = 1; i <= ntoy; ++i) {
    double x = gRandom->Gaus(xpeak, sigma);
    if (x >= xmin && x < xmax)
      hist->Fill(x);
  }
  return hist;
}

#pragma GCC diagnostic pop
