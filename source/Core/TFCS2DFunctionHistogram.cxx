// Copyright (c) 2024 CERN for the benefit of the FastCaloSim project

#include <algorithm>
#include <iostream>

#include "FastCaloSim/Core/TFCS2DFunctionHistogram.h"

#include "TH2.h"

//=============================================
//======= TFCS2DFunctionHistogram =========
//=============================================
using namespace FastCaloSim::Core;

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
          MSG_WARNING("bin content is negative in histogram "
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
    MSG_ERROR("histogram " << hist->GetName() << " : " << hist->GetTitle()
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
