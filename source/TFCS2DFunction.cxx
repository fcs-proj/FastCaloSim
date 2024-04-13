#pragma GCC diagnostic push 
#pragma GCC diagnostic ignored "-Wfloat-equal"

/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#include "FastCaloSim/TFCS2DFunction.h"
#include "FastCaloSim/TFCS2DFunctionHistogram.h"
#include "FastCaloSim/TFCS2DFunctionTemplateHistogram.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TRandom.h"
#include "TFile.h"

#include <iostream>

//=============================================
//======= TFCS2DFunction =========
//=============================================

void TFCS2DFunction::rnd_to_fct(float value[], const float rnd[]) const {
  rnd_to_fct(value[0], value[1], rnd[0], rnd[1]);
}

//================================================================================================================================

double TFCS2DFunction::CheckAndIntegrate2DHistogram(
    const TH2 *hist, std::vector<double> &integral_vec, int &first, int &last) {
  ISF_FCS::MLogging logger;
  Int_t nbinsx = hist->GetNbinsX();
  Int_t nbinsy = hist->GetNbinsY();
  Int_t nbins = nbinsx * nbinsy;

  float integral = 0.0;
  float hint = hist->Integral();
  integral_vec.resize(nbins);

  for (int ix = 1; ix <= nbinsx; ix++) {
    for (int iy = 1; iy <= nbinsy; iy++) {
      int globalbin = (ix - 1) * nbinsy + iy - 1;
      float binval = hist->GetBinContent(ix, iy);
      if (binval < 0) {
        // Can't work if a bin is negative, forcing bins to 0 in this case
        double fraction = binval / hint;
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
      integral_vec[globalbin] = integral;
    }
  }

  for (first = 0; first < nbins; first++)
    if (integral_vec[first] != 0)
      break;
  for (last = nbins - 1; last > 0; last--)
    if (integral_vec[last] != integral)
      break;
  last++;

  if (integral <= 0) {
    ATH_MSG_NOCLASS(logger, "Error: histogram "
                                << hist->GetName() << " : " << hist->GetTitle()
                                << " integral=" << integral << " is <=0");
  }

  return integral;
}

TH2 *create_random_TH2(int nbinsx = 64, int nbinsy = 64) {
  TH2 *hist = new TH2F("test2D", "test2D", nbinsx, 0, 1, nbinsy, 0, 1);
  hist->Sumw2();
  for (int ix = 1; ix <= nbinsx; ++ix) {
    for (int iy = 1; iy <= nbinsy; ++iy) {
      hist->SetBinContent(ix, iy,
                          (0.5 + gRandom->Rndm()) * (nbinsx + ix) *
                              (nbinsy * nbinsy / 2 + iy * iy));
      if (gRandom->Rndm() < 0.1)
        hist->SetBinContent(ix, iy, 0);
      hist->SetBinError(ix, iy, 0);
    }
  }
  return hist;
}

void TFCS2DFunction::unit_test(TH2 *hist,
                                                     TFCS2DFunction *rtof,
                                                     const char *outfilename,
                                                     int nrnd) {
  ISF_FCS::MLogging logger;
  if (hist == nullptr)
    hist = create_random_TH2();
  if (rtof == nullptr)
    rtof = new TFCS2DFunctionHistogram(hist);

  int nbinsx = hist->GetNbinsX();
  int nbinsy = hist->GetNbinsY();

  float value[2];
  float rnd[2];
  //cppcheck-suppress uninitvar
  for (rnd[0] = 0; rnd[0] < 0.9999; rnd[0] += 0.25) {
    for (rnd[1] = 0; rnd[1] < 0.9999; rnd[1] += 0.25) {
      rtof->rnd_to_fct(value, rnd);
      ATH_MSG_NOCLASS(logger, "rnd0=" << rnd[0] << " rnd1=" << rnd[1]
                                      << " -> x=" << value[0]
                                      << " y=" << value[1]);
    }
  }

  //  TH2F* hist_val=new
  //  TH2F("val2D","val2D",16,hist->GetXaxis()->GetXmin(),hist->GetXaxis()->GetXmax(),
  //                                          16,hist->GetYaxis()->GetXmin(),hist->GetYaxis()->GetXmax());
  TH2F *hist_val =
      (TH2F *)hist->Clone(TString(hist->GetName()) + "_" + rtof->ClassName());
  hist_val->Reset();

  float weight = hist->Integral() / nrnd;
  hist_val->Sumw2();
  for (int i = 0; i < nrnd; ++i) {
    rnd[0] = gRandom->Rndm();
    rnd[1] = gRandom->Rndm();
    rtof->rnd_to_fct(value, rnd);
    hist_val->Fill(value[0], value[1], weight);
  }
  hist_val->Add(hist, -1);

  TH1F *hist_pull =
      new TH1F(TString("pull_") + rtof->ClassName(),
               TString("pull for ") + rtof->ClassName(), 80, -4, 4);
  for (int ix = 1; ix <= nbinsx; ++ix) {
    for (int iy = 1; iy <= nbinsy; ++iy) {
      float val = hist_val->GetBinContent(ix, iy);
      float err = hist_val->GetBinError(ix, iy);
      if (err > 0)
        hist_pull->Fill(val / err);
      // ATH_MSG_NOCLASS(logger,"val="<<val<<" err="<<err);
    }
  }

  std::unique_ptr<TFile> outputfile(TFile::Open(outfilename, "UPDATE"));
  if (outputfile != nullptr) {
    hist->Write();
    hist_val->Write();
    hist_pull->Write();
    outputfile->ls();
  }

  new TCanvas(hist->GetName(), hist->GetTitle());
  hist->Draw("colz");

  new TCanvas(hist_val->GetName(), hist_val->GetTitle());
  hist_val->Draw("colz");

  new TCanvas(hist_pull->GetName(), hist_pull->GetTitle());
  hist_pull->Draw();
}

void TFCS2DFunction::unit_tests(TH2 *hist, const char *outfilename, int nrnd) {
  if (hist == nullptr)
    hist = create_random_TH2(16, 16);

  const int ntest = 4;
  TFCS2DFunction *tests[ntest];
  tests[0] = new TFCS2DFunctionHistogram(hist);
  tests[1] = new TFCS2DFunctionInt8Int8Int8Histogram(hist);
  tests[2] = new TFCS2DFunctionInt8Int8Int16Histogram(hist);
  tests[3] = new TFCS2DFunctionInt8Int8Int32Histogram(hist);

  for (int i = 0; i < ntest; ++i) {
    unit_test(hist, tests[i], outfilename, nrnd);
  }
}

#pragma GCC diagnostic pop