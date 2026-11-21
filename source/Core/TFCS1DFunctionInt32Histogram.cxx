// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>

#include "FastCaloSim/Core/TFCS1DFunctionInt32Histogram.h"

#include "FastCaloSim/Definitions/Index.h"
#include "TH1.h"

//=============================================
//======= TFCS1DFunctionInt32Histogram =========
//=============================================

const TFCS1DFunctionInt32Histogram::HistoContent_t
    TFCS1DFunctionInt32Histogram::s_MaxValue = UINT32_MAX;

void TFCS1DFunctionInt32Histogram::Initialize(const TH1* hist)
{
  const Int_t nbinsx = hist->GetNbinsX();
  const auto nbins = static_cast<std::size_t>(nbinsx);

  double integral = 0;
  m_HistoBorders.resize(nbins + 1);
  m_HistoContents.resize(nbins);
  std::vector<double> temp_HistoContents(nbins);
  std::size_t ibin = 0;
  for (int ix = 1; ix <= nbinsx; ix++) {
    double binval = hist->GetBinContent(ix);
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
    temp_HistoContents[ibin] = integral;
    ++ibin;
  }
  if (integral <= 0) {
    FCS_MSG_ERROR("histogram " << hist->GetName() << " : " << hist->GetTitle()
                               << " integral=" << integral << " is <=0");
    m_HistoBorders.resize(0);
    m_HistoContents.resize(0);
    return;
  }

  // The bin borders are stored in single precision
  for (std::size_t i = 0; i < nbins; i++) {
    m_HistoBorders[i] = static_cast<float>(
        hist->GetXaxis()->GetBinLowEdge(static_cast<int>(i) + 1));
  }
  m_HistoBorders[nbins] = static_cast<float>(hist->GetXaxis()->GetXmax());

  for (ibin = 0; ibin < nbins; ++ibin) {
    m_HistoContents[ibin] = static_cast<HistoContent_t>(
        s_MaxValue * (temp_HistoContents[ibin] / integral));
  }
}

auto TFCS1DFunctionInt32Histogram::rnd_to_fct(double rnd) const -> double
{
  if (m_HistoContents.empty()) {
    return 0;
  }
  const auto int_rnd = static_cast<HistoContent_t>(s_MaxValue * rnd);
  auto it =
      std::upper_bound(m_HistoContents.begin(), m_HistoContents.end(), int_rnd);
  auto ibin = FastCaloSim::index_of(m_HistoContents, it);
  if (ibin >= m_HistoContents.size()) {
    ibin = m_HistoContents.size() - 1;
  }

  HistoContent_t basecont = 0;
  if (ibin > 0) {
    basecont = m_HistoContents[ibin - 1];
  }

  const auto low_edge = static_cast<double>(m_HistoBorders[ibin]);
  const auto bin_width =
      static_cast<double>(m_HistoBorders[ibin + 1]) - low_edge;
  const HistoContent_t dcont = m_HistoContents[ibin] - basecont;
  if (dcont > 0) {
    return low_edge
        + bin_width * static_cast<double>(int_rnd - basecont)
        / static_cast<double>(dcont);
  }
  return low_edge + bin_width / 2;
}

auto TFCS1DFunctionInt32Histogram::operator==(const TFCS1DFunction& ref) const
    -> bool
{
  if (IsA() != ref.IsA()) {
    return false;
  }
  const auto& ref_typed = static_cast<const TFCS1DFunctionInt32Histogram&>(ref);

  if (m_HistoBorders != ref_typed.m_HistoBorders) {
    return false;
  }
  return m_HistoContents == ref_typed.m_HistoContents;
}
