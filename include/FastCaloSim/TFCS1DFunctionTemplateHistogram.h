#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-equal"
/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#ifndef ISF_FASTCALOSIMEVENT_TFCS1DFunctionTemplateHistogram_h
#  define ISF_FASTCALOSIMEVENT_TFCS1DFunctionTemplateHistogram_h

#  include <iostream>

#  include "FastCaloSim/TFCS1DFunction.h"
#  include "FastCaloSim/TFCS1DFunctionTemplateHelpers.h"
#  include "TH1.h"

template<typename Txvec, typename Ty, typename Trandom = float>
class TFCS1DFunctionTemplateHistogram : public TFCS1DFunction
{
public:
  typedef TFCS1DFunction_size_t size_t;
  typedef Trandom random_type;
  typedef Txvec xvec_type;
  typedef Ty y_value_type;

  TFCS1DFunctionTemplateHistogram(TH1* hist = nullptr, bool doprint = true)
  {
    if (hist)
      Initialize(hist, doprint);
  };
  ~TFCS1DFunctionTemplateHistogram() {};

  std::size_t MemorySizeArray() const
  {
    return m_HistoBorders.MemorySizeArray() + m_HistoContents.MemorySizeArray();
  };
  std::size_t MemorySize() const { return sizeof(*this) + MemorySizeArray(); };

  /// set number of bins
  void set_nbins(size_t nbins)
  {
    m_HistoBorders.set_nbins(nbins);
    m_HistoContents.set_nbins(nbins);
  };

  /// return number of bins
  inline size_t get_nbins() const { return m_HistoContents.get_nbins(); };

  /// Initialize from root histogram. Depending on the precision of the x- and
  /// y-axis, bins are merged if numerical identical
  void Initialize(const TH1* hist, bool doprint = true)
  {
    Int_t nbins = hist->GetNbinsX();

    std::vector<double> temp_HistoContents;
    int first, last;
    double integral =
        CheckAndIntegrate1DHistogram(hist, temp_HistoContents, first, last);
    if (integral <= 0) {
      set_nbins(0);
    } else {
      set_nbins(nbins);

      m_HistoBorders.SetMinMax(hist->GetXaxis()->GetBinLowEdge(first + 1),
                               hist->GetXaxis()->GetBinUpEdge(last + 1));
      Int_t ihist = 0;
      for (Int_t ibin = first; ibin <= last; ++ibin) {
        m_HistoBorders.SetBinLowEdge(ihist,
                                     hist->GetXaxis()->GetBinLowEdge(ibin + 1));
        if (ihist > 0)
          if (m_HistoBorders.GetBinLowEdge(ihist - 1)
              == m_HistoBorders.GetBinLowEdge(ihist))
          {
            if (doprint)
              ATH_MSG_INFO("Skip bin="
                           << ibin + 1
                           << " x=" << hist->GetXaxis()->GetBinLowEdge(ibin + 1)
                           << " fx=" << m_HistoBorders.GetBinLowEdge(ihist));
            --ihist;
            if (doprint)
              ATH_MSG_INFO("     bin="
                           << ibin
                           << " x=" << hist->GetXaxis()->GetBinLowEdge(ibin)
                           << " fx=" << m_HistoBorders.GetBinLowEdge(ihist));
          }
        m_HistoContents.set_fraction(ihist,
                                     temp_HistoContents[ibin] / integral);
        if (ihist > 0)
          if (m_HistoContents.get_fraction(ihist - 1)
              == m_HistoContents.get_fraction(ihist))
          {
            if (doprint)
              ATH_MSG_INFO("Skip fbin="
                           << ihist
                           << " fx=" << m_HistoBorders.GetBinLowEdge(ihist)
                           << " frac=" << m_HistoContents.get_fraction(ihist));
            --ihist;
            if (doprint)
              ATH_MSG_INFO("     fbin="
                           << ihist
                           << " fx=" << m_HistoBorders.GetBinLowEdge(ihist)
                           << " frac=" << m_HistoContents.get_fraction(ihist));
          }

        if (doprint)
          ATH_MSG_INFO("bin=" << ibin + 1 << " fbin=" << ihist << "/"
                              << m_HistoBorders.get_nbins() << " x=["
                              << hist->GetXaxis()->GetBinLowEdge(ibin + 1)
                              << "," << hist->GetXaxis()->GetBinUpEdge(ibin + 1)
                              << "] fx=" << m_HistoBorders.GetBinLowEdge(ihist)
                              << " int=" << temp_HistoContents[ibin] / integral
                              << " frac="
                              << m_HistoContents.get_fraction(ihist));

        ++ihist;
      }
      if (ihist != nbins) {
        set_nbins(ihist);
        m_HistoBorders.SetMinMax(hist->GetXaxis()->GetBinLowEdge(first + 1),
                                 hist->GetXaxis()->GetBinUpEdge(last + 1));
      }
    }
  }

  using TFCS1DFunction::rnd_to_fct;

  /// Function gets random number rnd in the range [0,1) as argument
  /// and returns function value according to a histogram distribution
  virtual double rnd_to_fct(double rnd) const
  {
    if (m_HistoContents.get_nbins() == 0)
      return 0;
    Trandom residual_rnd;
    size_t ibin = m_HistoContents.get_bin(rnd, residual_rnd);
    // ATH_MSG_INFO( fx="<<m_HistoBorders.GetBinLowEdge(ibin)<<"
    // frac="<<m_HistoContents.get_fraction(ibin)<<"
    // residual_rnd="<<residual_rnd);
    return m_HistoBorders.position(ibin, residual_rnd);
  }

  inline const Txvec& get_HistoBordersx() const { return m_HistoBorders; };
  inline Txvec& get_HistoBordersx() { return m_HistoBorders; };

  inline const TFCS1DFunction_HistogramContent<Ty, Trandom>& get_HistoContents()
      const
  {
    return m_HistoContents;
  };
  inline TFCS1DFunction_HistogramContent<Ty, Trandom>& get_HistoContents()
  {
    return m_HistoContents;
  };

protected:
  Txvec m_HistoBorders;
  TFCS1DFunction_HistogramContent<Ty, Trandom> m_HistoContents;

private:
  ClassDef(TFCS1DFunctionTemplateHistogram,
           1)  // TFCS1DFunctionTemplateHistogram
};

class TFCS1DFunctionInt8Int8Histogram
    : public TFCS1DFunctionTemplateHistogram<
          TFCS1DFunction_HistogramInt8BinEdges,
          uint8_t,
          float>
{
public:
  TFCS1DFunctionInt8Int8Histogram(TH1* h = nullptr)
      : TFCS1DFunctionTemplateHistogram<TFCS1DFunction_HistogramInt8BinEdges,
                                        uint8_t,
                                        float>(h) {};

  ClassDef(TFCS1DFunctionInt8Int8Histogram,
           1)  // TFCS1DFunctionInt8Int8Histogram
};

class TFCS1DFunctionInt8Int16Histogram
    : public TFCS1DFunctionTemplateHistogram<
          TFCS1DFunction_HistogramInt8BinEdges,
          uint16_t,
          float>
{
public:
  TFCS1DFunctionInt8Int16Histogram(TH1* h = nullptr)
      : TFCS1DFunctionTemplateHistogram<TFCS1DFunction_HistogramInt8BinEdges,
                                        uint16_t,
                                        float>(h) {};

  ClassDef(TFCS1DFunctionInt8Int16Histogram,
           1)  // TFCS1DFunctionInt8Int16Histogram
};

class TFCS1DFunctionInt8Int32Histogram
    : public TFCS1DFunctionTemplateHistogram<
          TFCS1DFunction_HistogramInt8BinEdges,
          uint32_t,
          float>
{
public:
  TFCS1DFunctionInt8Int32Histogram(TH1* h = nullptr)
      : TFCS1DFunctionTemplateHistogram<TFCS1DFunction_HistogramInt8BinEdges,
                                        uint32_t,
                                        float>(h) {};

  ClassDef(TFCS1DFunctionInt8Int32Histogram,
           1)  // TFCS1DFunctionInt8Int32Histogram
};

class TFCS1DFunctionInt16Int16Histogram
    : public TFCS1DFunctionTemplateHistogram<
          TFCS1DFunction_HistogramInt16BinEdges,
          uint16_t,
          float>
{
public:
  TFCS1DFunctionInt16Int16Histogram(TH1* h = nullptr)
      : TFCS1DFunctionTemplateHistogram<TFCS1DFunction_HistogramInt16BinEdges,
                                        uint16_t,
                                        float>(h) {};

  ClassDef(TFCS1DFunctionInt16Int16Histogram,
           1)  // TFCS1DFunctionInt16Int16Histogram
};

class TFCS1DFunctionInt16Int32Histogram
    : public TFCS1DFunctionTemplateHistogram<
          TFCS1DFunction_HistogramInt16BinEdges,
          uint32_t,
          float>
{
public:
  TFCS1DFunctionInt16Int32Histogram(TH1* h = nullptr)
      : TFCS1DFunctionTemplateHistogram<TFCS1DFunction_HistogramInt16BinEdges,
                                        uint32_t,
                                        float>(h) {};

  ClassDef(TFCS1DFunctionInt16Int32Histogram,
           1)  // TFCS1DFunctionInt16Int32Histogram
};

class TFCS1DFunctionInt32Int32Histogram
    : public TFCS1DFunctionTemplateHistogram<
          TFCS1DFunction_HistogramInt32BinEdges,
          uint32_t,
          float>
{
public:
  TFCS1DFunctionInt32Int32Histogram(TH1* h = nullptr)
      : TFCS1DFunctionTemplateHistogram<TFCS1DFunction_HistogramInt32BinEdges,
                                        uint32_t,
                                        float>(h) {};

  ClassDef(TFCS1DFunctionInt32Int32Histogram,
           1)  // TFCS1DFunctionInt32Int32Histogram
};

#endif

#pragma GCC diagnostic pop