// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#pragma once

#include <vector>

#include "FastCaloSim/Core/TFCS1DFunction.h"
#include "FastCaloSim/Core/TFCS1DFunctionTemplateHelpers.h"
#include "FastCaloSim/Definitions/FloatCompare.h"
#include "TH1.h"

template<typename Txvec, typename Ty, typename Trandom = float>
class TFCS1DFunctionTemplateHistogram : public TFCS1DFunction
{
public:
  using size_t = TFCS1DFunction_size_t;
  using random_type = Trandom;
  using xvec_type = Txvec;
  using y_value_type = Ty;

  TFCS1DFunctionTemplateHistogram(TH1* hist = nullptr, bool doprint = true)
  {
    if (hist) {
      Initialize(hist, doprint);
    }
  }
  ~TFCS1DFunctionTemplateHistogram() override = default;

  auto MemorySizeArray() const -> std::size_t
  {
    return m_HistoBorders.MemorySizeArray() + m_HistoContents.MemorySizeArray();
  }
  auto MemorySize() const -> std::size_t override
  {
    return sizeof(*this) + MemorySizeArray();
  }

  /// set number of bins
  void set_nbins(size_t nbins)
  {
    m_HistoBorders.set_nbins(nbins);
    m_HistoContents.set_nbins(nbins);
  }

  /// return number of bins
  auto get_nbins() const -> size_t { return m_HistoContents.get_nbins(); }

  /// Initialize from root histogram. Depending on the precision of the x- and
  /// y-axis, bins are merged if numerical identical
  void Initialize(const TH1* hist, bool doprint = true)
  {
    const auto nbins = static_cast<size_t>(hist->GetNbinsX());

    std::vector<double> temp_HistoContents;
    int first;
    int last;
    double const integral =
        CheckAndIntegrate1DHistogram(hist, temp_HistoContents, first, last);
    if (integral <= 0) {
      set_nbins(0);
      return;
    }
    set_nbins(nbins);

    // The bin edges are stored with the precision of Txvec
    using edge_type = typename Txvec::value_type;
    const TAxis* axis = hist->GetXaxis();
    const auto low_edge = [axis](int bin) -> auto
    { return static_cast<edge_type>(axis->GetBinLowEdge(bin)); };
    const auto up_edge = [axis](int bin) -> auto
    { return static_cast<edge_type>(axis->GetBinUpEdge(bin)); };

    m_HistoBorders.SetMinMax(low_edge(first + 1), up_edge(last + 1));
    size_t ihist = 0;
    for (int ibin = first; ibin <= last; ++ibin) {
      m_HistoBorders.SetBinLowEdge(ihist, low_edge(ibin + 1));
      // Merge bins that are numerically identical at the stored precision
      if (ihist > 0
          && FastCaloSim::exactly_equal(m_HistoBorders.GetBinLowEdge(ihist - 1),
                                        m_HistoBorders.GetBinLowEdge(ihist)))
      {
        if (doprint) {
          FCS_MSG_INFO("Skip bin="
                       << ibin + 1 << " x=" << axis->GetBinLowEdge(ibin + 1)
                       << " fx=" << m_HistoBorders.GetBinLowEdge(ihist));
        }
        --ihist;
        if (doprint) {
          FCS_MSG_INFO("     bin=" << ibin << " x=" << axis->GetBinLowEdge(ibin)
                                   << " fx="
                                   << m_HistoBorders.GetBinLowEdge(ihist));
        }
      }
      const double fraction =
          temp_HistoContents[static_cast<std::size_t>(ibin)] / integral;
      m_HistoContents.set_fraction(ihist, static_cast<Trandom>(fraction));
      if (ihist > 0
          && FastCaloSim::exactly_equal(m_HistoContents.get_fraction(ihist - 1),
                                        m_HistoContents.get_fraction(ihist)))
      {
        if (doprint) {
          FCS_MSG_INFO("Skip fbin="
                       << ihist << " fx=" << m_HistoBorders.GetBinLowEdge(ihist)
                       << " frac=" << m_HistoContents.get_fraction(ihist));
        }
        --ihist;
        if (doprint) {
          FCS_MSG_INFO("     fbin="
                       << ihist << " fx=" << m_HistoBorders.GetBinLowEdge(ihist)
                       << " frac=" << m_HistoContents.get_fraction(ihist));
        }
      }

      if (doprint) {
        FCS_MSG_INFO("bin=" << ibin + 1 << " fbin=" << ihist << "/"
                            << m_HistoBorders.get_nbins() << " x=["
                            << axis->GetBinLowEdge(ibin + 1) << ","
                            << axis->GetBinUpEdge(ibin + 1)
                            << "] fx=" << m_HistoBorders.GetBinLowEdge(ihist)
                            << " int=" << fraction
                            << " frac=" << m_HistoContents.get_fraction(ihist));
      }

      ++ihist;
    }
    if (ihist != nbins) {
      set_nbins(ihist);
      m_HistoBorders.SetMinMax(low_edge(first + 1), up_edge(last + 1));
    }
  }

  using TFCS1DFunction::rnd_to_fct;

  /// Function gets random number rnd in the range [0,1) as argument
  /// and returns function value according to a histogram distribution
  auto rnd_to_fct(double rnd) const -> double override
  {
    if (m_HistoContents.get_nbins() == 0) {
      return 0;
    }
    Trandom residual_rnd;
    size_t const ibin =
        m_HistoContents.get_bin(static_cast<Trandom>(rnd), residual_rnd);
    return static_cast<double>(m_HistoBorders.position(ibin, residual_rnd));
  }

  auto get_HistoBordersx() const -> const Txvec& { return m_HistoBorders; }
  auto get_HistoBordersx() -> Txvec& { return m_HistoBorders; }

  auto get_HistoContents() const
      -> const TFCS1DFunction_HistogramContent<Ty, Trandom>&
  {
    return m_HistoContents;
  }
  auto get_HistoContents() -> TFCS1DFunction_HistogramContent<Ty, Trandom>&
  {
    return m_HistoContents;
  }

protected:
  Txvec m_HistoBorders;
  TFCS1DFunction_HistogramContent<Ty, Trandom> m_HistoContents;

private:
  ClassDefOverride(TFCS1DFunctionTemplateHistogram,
                   1)  // TFCS1DFunctionTemplateHistogram
};

class TFCS1DFunctionInt8Int8Histogram
    : public TFCS1DFunctionTemplateHistogram<
          TFCS1DFunction_HistogramInt8BinEdges,
          uint8_t>
{
public:
  using TFCS1DFunctionTemplateHistogram::TFCS1DFunctionTemplateHistogram;

  ClassDefOverride(TFCS1DFunctionInt8Int8Histogram,
                   1)  // TFCS1DFunctionInt8Int8Histogram
};

class TFCS1DFunctionInt8Int16Histogram
    : public TFCS1DFunctionTemplateHistogram<
          TFCS1DFunction_HistogramInt8BinEdges,
          uint16_t>
{
public:
  using TFCS1DFunctionTemplateHistogram::TFCS1DFunctionTemplateHistogram;

  ClassDefOverride(TFCS1DFunctionInt8Int16Histogram,
                   1)  // TFCS1DFunctionInt8Int16Histogram
};

class TFCS1DFunctionInt8Int32Histogram
    : public TFCS1DFunctionTemplateHistogram<
          TFCS1DFunction_HistogramInt8BinEdges,
          uint32_t>
{
public:
  using TFCS1DFunctionTemplateHistogram::TFCS1DFunctionTemplateHistogram;

  ClassDefOverride(TFCS1DFunctionInt8Int32Histogram,
                   1)  // TFCS1DFunctionInt8Int32Histogram
};

class TFCS1DFunctionInt32Int32Histogram
    : public TFCS1DFunctionTemplateHistogram<
          TFCS1DFunction_HistogramInt32BinEdges,
          uint32_t>
{
public:
  using TFCS1DFunctionTemplateHistogram::TFCS1DFunctionTemplateHistogram;

  ClassDefOverride(TFCS1DFunctionInt32Int32Histogram,
                   1)  // TFCS1DFunctionInt32Int32Histogram
};
