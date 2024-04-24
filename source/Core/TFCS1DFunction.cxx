#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-equal"
/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#include <iostream>
#include <string>

#include "FastCaloSim/Core/TFCS1DFunction.h"

#include "TCanvas.h"
#include "TH1.h"
#include "TRandom.h"

//=============================================
//======= TFCS1DFunction =========
//=============================================

void TFCS1DFunction::rnd_to_fct(float value[], const float rnd[]) const
{
  value[0] = rnd_to_fct(rnd[0]);
}

double TFCS1DFunction::get_maxdev(TH1* h_input1, TH1* h_approx1)
{
  TH1D* h_input = (TH1D*)h_input1->Clone("h_input");
  TH1D* h_approx = (TH1D*)h_approx1->Clone("h_approx");

  double maxdev = 0.0;

  // normalize the histos to the same area:
  double integral_input = h_input->Integral();
  double integral_approx = 0.0;
  for (int b = 1; b <= h_input->GetNbinsX(); b++)
    integral_approx +=
        h_approx->GetBinContent(h_approx->FindBin(h_input->GetBinCenter(b)));
  h_approx->Scale(integral_input / integral_approx);

  double ymax = h_approx->GetBinContent(h_approx->GetNbinsX())
      - h_approx->GetBinContent(h_approx->GetMinimumBin());
  for (int i = 1; i <= h_input->GetNbinsX(); i++) {
    double val = fabs(h_approx->GetBinContent(
                          h_approx->FindBin(h_input->GetBinCenter(i)))
                      - h_input->GetBinContent(i))
        / ymax;
    if (val > maxdev)
      maxdev = val;
  }

  delete h_input;
  delete h_approx;

  return maxdev * 100.0;
}

double TFCS1DFunction::CheckAndIntegrate1DHistogram(
    const TH1* hist, std::vector<double>& integral_vec, int& first, int& last)
{
  ISF_FCS::MLogging logger;
  Int_t nbins = hist->GetNbinsX();

  float integral = 0;
  integral_vec.resize(nbins);
  for (int ix = 1; ix <= nbins; ix++) {
    float binval = hist->GetBinContent(ix);
    if (binval < 0) {
      // Can't work if a bin is negative, forcing bins to 0 in this case
      double fraction = binval / hist->Integral();
      if (TMath::Abs(fraction) > 1e-5) {
        ATH_MSG_NOCLASS(logger,
                        "Warning: bin content is negative in histogram "
                            << hist->GetName() << " : " << hist->GetTitle()
                            << " binval=" << binval << " " << fraction * 100
                            << "% of integral=" << hist->Integral()
                            << ". Forcing bin to 0.");
      }
      binval = 0;
    }
    integral += binval;
    integral_vec[ix - 1] = integral;
  }

  for (first = 0; first < nbins; first++)
    if (integral_vec[first] != 0)
      break;
  for (last = nbins - 1; last > 0; last--)
    if (integral_vec[last] != integral)
      break;
  last++;

  if (integral <= 0) {
    ATH_MSG_NOCLASS(logger,
                    "Error: histogram "
                        << hist->GetName() << " : " << hist->GetTitle()
                        << " integral=" << integral << " is <=0");
  }
  return integral;
}

TH1* TFCS1DFunction::generate_histogram_random_slope(int nbinsx,
                                                     double xmin,
                                                     double xmax,
                                                     double zerothreshold)
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

TH1* TFCS1DFunction::generate_histogram_random_gauss(
    int nbinsx, int ntoy, double xmin, double xmax, double xpeak, double sigma)
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

void TFCS1DFunction::unit_test(TH1* hist,
                               TFCS1DFunction* rtof,
                               int nrnd,
                               TH1* histfine)
{
  ISF_FCS::MLogging logger;
  ATH_MSG_NOCLASS(logger,
                  "========= " << hist->GetName() << " funcsize="
                               << rtof->MemorySize() << " ========");
  int nbinsx = hist->GetNbinsX();
  double integral = hist->Integral();

  float value[2];
  float rnd[2];
  // cppcheck-suppress uninitvar
  for (rnd[0] = 0; rnd[0] < 0.9999; rnd[0] += 0.25) {
    rtof->rnd_to_fct(value, rnd);
    ATH_MSG_NOCLASS(logger, "rnd0=" << rnd[0] << " -> x=" << value[0]);
  }

  TH1* hist_val = nullptr;
  if (histfine)
    hist_val = (TH1*)histfine->Clone(TString(hist->GetName()) + "hist_val");
  else
    hist_val = (TH1*)hist->Clone(TString(hist->GetName()) + "hist_val");
  double weightfine = hist_val->Integral() / nrnd;
  hist_val->SetTitle("toy simulation");
  hist_val->Reset();
  hist_val->SetLineColor(2);
  hist_val->Sumw2();

  TH1* hist_diff = (TH1*)hist->Clone(TString(hist->GetName()) + "_difference");
  hist_diff->SetTitle("cut efficiency difference");
  hist_diff->Reset();
  hist_diff->Sumw2();

  double weight = integral / nrnd;
  for (int i = 0; i < nrnd; ++i) {
    rnd[0] = gRandom->Rndm();
    rtof->rnd_to_fct(value, rnd);
    hist_val->Fill(value[0], weightfine);
    hist_diff->Fill(value[0], weight);
  }
  hist_diff->Add(hist, -1);
  hist_diff->Scale(1.0 / integral);

  TH1F* hist_pull =
      new TH1F(TString(hist->GetName()) + "_pull", "pull", 200, -10, 10);
  for (int ix = 1; ix <= nbinsx; ++ix) {
    float val = hist_diff->GetBinContent(ix);
    float err = hist_diff->GetBinError(ix);
    if (err > 0)
      hist_pull->Fill(val / err);
  }

  TCanvas* c = new TCanvas(hist->GetName(), hist->GetName());
  c->Divide(2, 2);

  c->cd(1);
  if (histfine) {
    histfine->SetLineColor(kGray);
    histfine->DrawClone("hist");
    hist->DrawClone("same");
  } else {
    hist->DrawClone();
  }
  hist_val->DrawClone("sameshist");

  c->cd(2);
  if (histfine) {
    histfine->SetLineColor(kGray);
    histfine->DrawClone("hist");
    hist->DrawClone("same");
  } else {
    hist->DrawClone();
  }
  hist_val->DrawClone("sameshist");
  gPad->SetLogy();

  c->cd(3);
  hist_diff->Draw();

  c->cd(4);
  hist_pull->Draw();

  c->SaveAs(".png");
}

#pragma GCC diagnostic pop
