/*
  Copyright (C) 2002-2023 CERN for the benefit of the ATLAS collaboration
*/

#include <algorithm>
#include <iostream>

#include "FastCaloSim/TFCS2DFunctionHistogram.h"

#include "TCanvas.h"
#include "TFile.h"
#include "TH2F.h"
#include "TMath.h"
#include "TRandom.h"

//=============================================
//======= TFCS2DFunctionHistogram =========
//=============================================

void TFCS2DFunctionHistogram::Initialize(TH2* hist)
{
  Int_t nbinsx = hist->GetNbinsX();
  Int_t nbinsy = hist->GetNbinsY();
  Int_t nbins = nbinsx * nbinsy;

  float integral = 0;
  m_HistoBorders.resize(nbinsx + 1);
  m_HistoBordersy.resize(nbinsy + 1);
  m_HistoContents.resize(nbins);
  int ibin = 0;
  for (int iy = 1; iy <= nbinsy; iy++) {
    for (int ix = 1; ix <= nbinsx; ix++) {
      float binval = hist->GetBinContent(ix, iy);
      if (binval < 0) {
        // Can't work if a bin is negative, forcing bins to 0 in this case
        double fraction = binval / hist->Integral();
        if (TMath::Abs(fraction) > 1e-5) {
          ATH_MSG_WARNING("bin content is negative in histogram "
                          << hist->GetName() << " : " << hist->GetTitle()
                          << " binval=" << binval << " " << fraction * 100
                          << "% of integral=" << hist->Integral()
                          << ". Forcing bin to 0.");
        }
        binval = 0;
      }
      integral += binval;
      m_HistoContents[ibin] = integral;
      ++ibin;
    }
  }
  if (integral <= 0) {
    ATH_MSG_ERROR("histogram " << hist->GetName() << " : " << hist->GetTitle()
                               << " integral=" << integral << " is <=0");
    m_HistoBorders.resize(0);
    m_HistoBordersy.resize(0);
    m_HistoContents.resize(0);
    return;
  }

  for (int ix = 1; ix <= nbinsx; ix++)
    m_HistoBorders[ix - 1] = hist->GetXaxis()->GetBinLowEdge(ix);
  m_HistoBorders[nbinsx] = hist->GetXaxis()->GetXmax();

  for (int iy = 1; iy <= nbinsy; iy++)
    m_HistoBordersy[iy - 1] = hist->GetYaxis()->GetBinLowEdge(iy);
  m_HistoBordersy[nbinsy] = hist->GetYaxis()->GetXmax();

  for (ibin = 0; ibin < nbins; ++ibin)
    m_HistoContents[ibin] /= integral;
}

void TFCS2DFunctionHistogram::rnd_to_fct(float& valuex,
                                         float& valuey,
                                         float rnd0,
                                         float rnd1) const
{
  if (m_HistoContents.empty()) {
    valuex = 0;
    valuey = 0;
    return;
  }
  auto it =
      std::upper_bound(m_HistoContents.begin(), m_HistoContents.end(), rnd0);
  int ibin = std::distance(m_HistoContents.begin(), it);
  if (ibin >= (int)m_HistoContents.size())
    ibin = m_HistoContents.size() - 1;
  Int_t nbinsx = m_HistoBorders.size() - 1;
  Int_t biny = ibin / nbinsx;
  Int_t binx = ibin - nbinsx * biny;

  float basecont = 0;
  if (ibin > 0)
    basecont = m_HistoContents[ibin - 1];

  float dcont = m_HistoContents[ibin] - basecont;
  if (dcont > 0) {
    valuex = m_HistoBorders[binx]
        + (m_HistoBorders[binx + 1] - m_HistoBorders[binx]) * (rnd0 - basecont)
            / dcont;
  } else {
    valuex = m_HistoBorders[binx]
        + (m_HistoBorders[binx + 1] - m_HistoBorders[binx]) / 2;
  }
  valuey = m_HistoBordersy[biny]
      + (m_HistoBordersy[biny + 1] - m_HistoBordersy[biny]) * rnd1;
}

void TFCS2DFunctionHistogram::unit_test(TH2* hist)
{
  ISF_FCS::MLogging logger;
  int nbinsx;
  int nbinsy;
  if (hist == nullptr) {
    //    hist=new TH2F("test2D","test2D",5,0,5,5,0,10);
    nbinsx = 64;
    nbinsy = 64;
    hist = new TH2F("test2D", "test2D", nbinsx, 0, 1, nbinsy, 0, 1);
    hist->Sumw2();
    for (int ix = 1; ix <= nbinsx; ++ix) {
      for (int iy = 1; iy <= nbinsy; ++iy) {
        hist->SetBinContent(ix,
                            iy,
                            (0.5 + gRandom->Rndm()) * (nbinsx + ix)
                                * (nbinsy * nbinsy / 2 + iy * iy));
        if (gRandom->Rndm() < 0.1)
          hist->SetBinContent(ix, iy, 0);
        hist->SetBinError(ix, iy, 0);
      }
    }
  }
  TFCS2DFunctionHistogram rtof(hist);
  nbinsx = hist->GetNbinsX();
  nbinsy = hist->GetNbinsY();

  float value[2];
  float rnd[2];
  // cppcheck-suppress uninitvar
  for (rnd[0] = 0; rnd[0] < 0.9999; rnd[0] += 0.25) {
    for (rnd[1] = 0; rnd[1] < 0.9999; rnd[1] += 0.25) {
      rtof.rnd_to_fct(value, rnd);
      ATH_MSG_NOCLASS(logger,
                      "rnd0=" << rnd[0] << " rnd1=" << rnd[1]
                              << " -> x=" << value[0] << " y=" << value[1]);
    }
  }

  //  TH2F* hist_val=new
  //  TH2F("val2D","val2D",16,hist->GetXaxis()->GetXmin(),hist->GetXaxis()->GetXmax(),
  //                                          16,hist->GetYaxis()->GetXmin(),hist->GetYaxis()->GetXmax());
  TH2F* hist_val = (TH2F*)hist->Clone("hist_val");
  hist_val->Reset();
  int nrnd = 100000000;
  float weight = hist->Integral() / nrnd;
  hist_val->Sumw2();
  for (int i = 0; i < nrnd; ++i) {
    rnd[0] = gRandom->Rndm();
    rnd[1] = gRandom->Rndm();
    rtof.rnd_to_fct(value, rnd);
    hist_val->Fill(value[0], value[1], weight);
  }
  hist_val->Add(hist, -1);

  TH1F* hist_pull = new TH1F("pull", "pull", 80, -4, 4);
  for (int ix = 1; ix <= nbinsx; ++ix) {
    for (int iy = 1; iy <= nbinsy; ++iy) {
      float val = hist_val->GetBinContent(ix, iy);
      float err = hist_val->GetBinError(ix, iy);
      if (err > 0)
        hist_pull->Fill(val / err);
      ATH_MSG_NOCLASS(logger, "val=" << val << " err=" << err);
    }
  }

  std::unique_ptr<TFile> outputfile(
      TFile::Open("TFCS2DFunctionHistogram_unit_test.root", "RECREATE"));
  if (outputfile != nullptr) {
    hist->Write();
    hist_val->Write();
    hist_pull->Write();
    outputfile->ls();
    outputfile->Close();
  }

  new TCanvas("input", "Input");
  hist->Draw("colz");

  new TCanvas("validation", "Validation");
  hist_val->Draw("colz");

  new TCanvas("pull", "Pull");
  hist_pull->Draw();
}
