// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#include <algorithm>
#include <cmath>
#include <iostream>

#include "FastCaloSim/Core/TFCS2DFunctionHistogram.h"

#include "FastCaloSim/Definitions/Index.h"
#include "TH2.h"

//=============================================
//======= TFCS2DFunctionHistogram =========
//=============================================

void TFCS2DFunctionHistogram::Initialize(TH2* hist)
{
  Int_t const nbinsx = hist->GetNbinsX();
  Int_t const nbinsy = hist->GetNbinsY();
  Int_t const nbins = nbinsx * nbinsy;

  // Accumulated in double precision, stored in single precision
  double integral = 0;
  const auto size_x = static_cast<std::size_t>(nbinsx);
  const auto size_y = static_cast<std::size_t>(nbinsy);
  m_HistoBorders.resize(size_x + 1);
  m_HistoBordersy.resize(size_y + 1);
  m_HistoContents.resize(static_cast<std::size_t>(nbins));
  std::size_t ibin = 0;
  for (int iy = 1; iy <= nbinsy; iy++) {
    for (int ix = 1; ix <= nbinsx; ix++) {
      double binval = hist->GetBinContent(ix, iy);
      if (binval < 0) {
        // Can't work if a bin is negative, forcing bins to 0 in this case
        double const fraction = binval / hist->Integral();
        if (std::abs(fraction) > 1e-5) {
          FCS_MSG_WARNING("bin content is negative in histogram "
                          << hist->GetName() << " : " << hist->GetTitle()
                          << " binval=" << binval << " " << fraction * 100
                          << "% of integral=" << hist->Integral()
                          << ". Forcing bin to 0.");
        }
        binval = 0;
      }
      integral += binval;
      m_HistoContents[ibin] = static_cast<float>(integral);
      ++ibin;
    }
  }
  if (integral <= 0) {
    FCS_MSG_ERROR("histogram " << hist->GetName() << " : " << hist->GetTitle()
                               << " integral=" << integral << " is <=0");
    m_HistoBorders.resize(0);
    m_HistoBordersy.resize(0);
    m_HistoContents.resize(0);
    return;
  }

  for (std::size_t ix = 0; ix < size_x; ix++) {
    m_HistoBorders[ix] = static_cast<float>(
        hist->GetXaxis()->GetBinLowEdge(static_cast<int>(ix) + 1));
  }
  m_HistoBorders[size_x] = static_cast<float>(hist->GetXaxis()->GetXmax());

  for (std::size_t iy = 0; iy < size_y; iy++) {
    m_HistoBordersy[iy] = static_cast<float>(
        hist->GetYaxis()->GetBinLowEdge(static_cast<int>(iy) + 1));
  }
  m_HistoBordersy[size_y] = static_cast<float>(hist->GetYaxis()->GetXmax());

  for (float& content : m_HistoContents) {
    content = static_cast<float>(static_cast<double>(content) / integral);
  }
}

void TFCS2DFunctionHistogram::rnd_to_fct(double& valuex,
                                         double& valuey,
                                         double rnd0,
                                         double rnd1) const
{
  if (m_HistoContents.empty()) {
    valuex = 0;
    valuey = 0;
    return;
  }
  // The bin borders and contents are stored in single precision
  const auto content = [this](std::size_t bin) -> double
  { return static_cast<double>(m_HistoContents[bin]); };
  const auto border_x = [this](std::size_t bin) -> double
  { return static_cast<double>(m_HistoBorders[bin]); };
  const auto border_y = [this](std::size_t bin) -> double
  { return static_cast<double>(m_HistoBordersy[bin]); };

  const auto it = std::ranges::upper_bound(m_HistoContents,
                                           rnd0,
                                           std::ranges::less {},
                                           [](float c) -> double
                                           { return static_cast<double>(c); });
  std::size_t ibin = FastCaloSim::index_of(m_HistoContents, it);
  if (ibin >= m_HistoContents.size()) {
    ibin = m_HistoContents.size() - 1;
  }
  const std::size_t nbinsx = m_HistoBorders.size() - 1;
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
