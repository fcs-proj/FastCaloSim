/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#include <algorithm>
#include <iostream>

#include "FastCaloSim/Core/TFCS1DFunctionInt16Histogram.h"

#include "TCanvas.h"
#include "TFile.h"
#include "TH2F.h"
#include "TMath.h"
#include "TRandom.h"

//=============================================
//======= TFCS1DFunctionInt16Histogram =========
//=============================================

const TFCS1DFunctionInt16Histogram::HistoContent_t
    TFCS1DFunctionInt16Histogram::s_MaxValue = UINT16_MAX;

void TFCS1DFunctionInt16Histogram::Initialize(const TH1* hist)
{
  Int_t nbinsx = hist->GetNbinsX();
  Int_t nbins = nbinsx;

  float integral = 0;
  m_HistoBorders.resize(nbinsx + 1);
  m_HistoContents.resize(nbins);
  std::vector<double> temp_HistoContents(nbins);
  int ibin = 0;
  for (int ix = 1; ix <= nbinsx; ix++) {
    float binval = hist->GetBinContent(ix);
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
    temp_HistoContents[ibin] = integral;
    ++ibin;
  }
  if (integral <= 0) {
    ATH_MSG_ERROR("histogram " << hist->GetName() << " : " << hist->GetTitle()
                               << " integral=" << integral << " is <=0");
    m_HistoBorders.resize(0);
    m_HistoContents.resize(0);
    return;
  }

  for (int ix = 1; ix <= nbinsx; ix++)
    m_HistoBorders[ix - 1] = hist->GetXaxis()->GetBinLowEdge(ix);
  m_HistoBorders[nbinsx] = hist->GetXaxis()->GetXmax();

  for (ibin = 0; ibin < nbins; ++ibin) {
    m_HistoContents[ibin] = s_MaxValue * (temp_HistoContents[ibin] / integral);
    // ATH_MSG_INFO("bin="<<ibin<<" val="<<m_HistoContents[ibin]);
  }
}

double TFCS1DFunctionInt16Histogram::rnd_to_fct(double rnd) const
{
  if (m_HistoContents.empty()) {
    return 0;
  }
  HistoContent_t int_rnd = s_MaxValue * rnd;
  auto it =
      std::upper_bound(m_HistoContents.begin(), m_HistoContents.end(), int_rnd);
  int ibin = std::distance(m_HistoContents.begin(), it);
  if (ibin >= (int)m_HistoContents.size())
    ibin = m_HistoContents.size() - 1;
  Int_t binx = ibin;

  HistoContent_t basecont = 0;
  if (ibin > 0)
    basecont = m_HistoContents[ibin - 1];

  HistoContent_t dcont = m_HistoContents[ibin] - basecont;
  if (dcont > 0) {
    return m_HistoBorders[binx]
        + ((m_HistoBorders[binx + 1] - m_HistoBorders[binx])
           * (int_rnd - basecont))
        / dcont;
  } else {
    return m_HistoBorders[binx]
        + (m_HistoBorders[binx + 1] - m_HistoBorders[binx]) / 2;
  }
}

void TFCS1DFunctionInt16Histogram::unit_test(TH1* hist)
{
  ISF_FCS::MLogging logger;
  int nbinsx;
  if (hist == nullptr) {
    nbinsx = 50;
    hist = new TH1D("test1D", "test1D", nbinsx, 0, 1);
    hist->Sumw2();
    for (int ix = 1; ix <= nbinsx; ++ix) {
      double val = (0.5 + gRandom->Rndm()) * (nbinsx + ix);
      if (gRandom->Rndm() < 0.1)
        val = 0;
      hist->SetBinContent(ix, val);
      hist->SetBinError(ix, 0);
    }
  }
  TFCS1DFunctionInt16Histogram rtof(hist);
  nbinsx = hist->GetNbinsX();

  float value[2];
  float rnd[2];
  // cppcheck-suppress uninitvar
  for (rnd[0] = 0; rnd[0] < 0.9999; rnd[0] += 0.25) {
    rtof.rnd_to_fct(value, rnd);
    ATH_MSG_NOCLASS(logger, "rnd0=" << rnd[0] << " -> x=" << value[0]);
  }

  TH1* hist_val = (TH1*)hist->Clone("hist_val");
  hist_val->SetTitle("difference");
  hist_val->Reset();
  int nrnd = 100000;
  double weight = hist->Integral() / nrnd;
  hist_val->Sumw2();
  for (int i = 0; i < nrnd; ++i) {
    rnd[0] = gRandom->Rndm();
    rtof.rnd_to_fct(value, rnd);
    hist_val->Fill(value[0], weight);
  }
  hist_val->Add(hist, -1);

  TH1F* hist_pull = new TH1F("pull", "pull", 200, -10, 10);
  for (int ix = 1; ix <= nbinsx; ++ix) {
    float val = hist_val->GetBinContent(ix);
    float err = hist_val->GetBinError(ix);
    if (err > 0)
      hist_pull->Fill(val / err);
    ATH_MSG_NOCLASS(logger, "val=" << val << " err=" << err);
  }

  new TCanvas("input", "Input");
  hist->Draw();

  new TCanvas("validation", "Validation");
  hist_val->Draw();

  new TCanvas("pull", "Pull");
  hist_pull->Draw();
}
