/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#ifndef ISF_FASTCALOSIMEVENT_TFCS2DFunctionTemplateHistogram_h
#define ISF_FASTCALOSIMEVENT_TFCS2DFunctionTemplateHistogram_h

#include <iostream>

#include "FastCaloSim/Core/TFCS1DFunctionTemplateHelpers.h"
#include "TBuffer.h"
#include "TFCS2DFunction.h"
#include "TH1.h"
#include "TH2.h"

template<typename Txvec, typename Tyvec, typename Tz, typename Trandom = float>
class TFCS2DFunctionTemplateHistogram : public TFCS2DFunction
{
public:
  typedef TFCS1DFunction_size_t size_t;
  typedef Trandom random_type;
  typedef Txvec xvec_type;
  typedef Tyvec yvec_type;
  typedef Tz z_value_type;

  TFCS2DFunctionTemplateHistogram(TH2* hist = nullptr)
  {
    if (hist)
      Initialize(hist);
  };
  ~TFCS2DFunctionTemplateHistogram() {};

  std::size_t MemorySizeArray() const
  {
    return m_HistoBordersx.MemorySizeArray() + m_HistoBordersy.MemorySizeArray()
        + m_HistoContents.MemorySizeArray();
  };
  std::size_t MemorySize() const { return sizeof(*this) + MemorySizeArray(); };

  /// set number of bins
  void set_nbins(size_t nbinsx, size_t nbinsy)
  {
    m_HistoBordersx.set_nbins(nbinsx);
    m_HistoBordersy.set_nbins(nbinsy);
    m_HistoContents.set_nbins(nbinsx * nbinsy);
  };

  /// return number of bins
  inline size_t get_nbins() const { return m_HistoContents.get_nbins(); };
  // inline size_t get_nbinsx() const {return m_HistoBordersx.get_nbins();};
  // inline size_t get_nbinsy() const {return m_HistoBordersy.get_nbins();};

  /// Initialize from root histogram. Depending on the precision of the x- and
  /// y-axis, bins are merged if numerical identical
  void Initialize(const TH2* hist)
  {
    Int_t nbinsx = hist->GetNbinsX();
    Int_t nbinsy = hist->GetNbinsY();
    std::vector<double> temp_HistoContents;
    int first, last;
    double integral =
        CheckAndIntegrate2DHistogram(hist, temp_HistoContents, first, last);
    printf("initializing 2D template...");
    if (integral <= 0) {
      set_nbins(0, 0);
    } else {
      set_nbins(nbinsx, nbinsy);

      // m_HistoBordersx.SetMinMax(hist->GetXaxis()->GetBinLowEdge(first+1),hist->GetXaxis()->GetBinUpEdge(last+1));
      m_HistoBordersx.SetMinMax(hist->GetXaxis()->GetBinLowEdge(1),
                                hist->GetXaxis()->GetBinUpEdge(nbinsx));
      m_HistoBordersy.SetMinMax(hist->GetYaxis()->GetBinLowEdge(1),
                                hist->GetYaxis()->GetBinUpEdge(nbinsy));
      // Int_t ihist=0;
      for (Int_t ibinx = 1; ibinx <= hist->GetNbinsX(); ++ibinx) {
        m_HistoBordersx.SetBinLowEdge(
            ibinx, hist->GetXaxis()->GetBinLowEdge(ibinx + 1));
        //          if(ihist>0)
        //          if(m_HistoBorders.GetBinLowEdge(ihist-1)==m_HistoBorders.GetBinLowEdge(ihist))
        //          { ATH_MSG_INFO( bin="<<ibin+1<<"
        //          x="<<hist->GetXaxis()->GetBinLowEdge(ibin+1)<<"
        //          fx="<<m_HistoBorders.GetBinLowEdge(ihist);
        //          --ihist;
        //          ATH_MSG_INFO("     bin="<<ibin  <<"
        //          x="<<hist->GetXaxis()->GetBinLowEdge(ibin  )<<"
        //          fx="<<m_HistoBorders.GetBinLowEdge(ihist));
        //          }
      }  // for
      for (Int_t ibiny = 1; ibiny <= hist->GetNbinsY(); ++ibiny) {
        m_HistoBordersy.SetBinLowEdge(
            ibiny, hist->GetYaxis()->GetBinLowEdge(ibiny + 1));
        //          if(ihist>0)
        //          if(m_HistoBorders.GetBinLowEdge(ihist-1)==m_HistoBorders.GetBinLowEdge(ihist))
        //          { ATH_MSG_INFO("Skip bin="<<ibin+1<<"
        //          x="<<hist->GetXaxis()->GetBinLowEdge(ibin+1)<<"
        //          fx="<<m_HistoBorders.GetBinLowEdge(ihist));
        //          --ihist;
        //          ATH_MSG_INFO("     bin="<<ibin  <<"
        //          x="<<hist->GetXaxis()->GetBinLowEdge(ibin  )<<"
        //          fx="<<m_HistoBorders.GetBinLowEdge(ihist));
        //          }
      }  // for

      for (Int_t ibinx = 1; ibinx <= hist->GetNbinsX(); ++ibinx)
        for (Int_t ibiny = 1; ibiny <= hist->GetNbinsY(); ++ibiny) {
          int globalbin = (ibinx - 1) * hist->GetNbinsY() + ibiny - 1;
          m_HistoContents.set_fraction(
              globalbin, temp_HistoContents[globalbin] / integral);
          //          if(ihist>0)
          //          if(m_HistoContents.get_fraction(ihist-1)==m_HistoContents.get_fraction(ihist))
          //          {
          //            ATH_MSG_INFO("Skip fbin="<<ihist<<"
          //            fx="<<m_HistoBorders.GetBinLowEdge(ihist)<<"
          //            frac="<<m_HistoContents.get_fraction(ihist));
          //            --ihist;
          //            ATH_MSG_INFO("     fbin="<<ihist<<"
          //            fx="<<m_HistoBorders.GetBinLowEdge(ihist)<<"
          //            frac="<<m_HistoContents.get_fraction(ihist));
          //          }

          //          ATH_MSG_INFO("bin="<<ibin+1<<"
          //          fbin="<<ihist<<"/"<<m_HistoBorders.get_nbins()<<"
          //          x=["<<hist->GetXaxis()->GetBinLowEdge(ibin+1)<<","<<hist->GetXaxis()->GetBinUpEdge(ibin+1)<<"]
          //          fx="<<m_HistoBorders.GetBinLowEdge(ihist)<<"
          //          int="<<temp_HistoContents[ibin]/integral<<"
          //          frac="<<m_HistoContents.get_fraction(ihist));

          //          ++ihist;
        }
      //        if(ihist!=nbins) {
      //          set_nbins(ihist);
      //          m_HistoBorders.SetMinMax(hist->GetXaxis()->GetBinLowEdge(first+1),hist->GetXaxis()->GetBinUpEdge(last+1));
      //        }
    }  // for ibinx
  }  // for ibiny

  /*
          for( Int_t ibiny=1; ibiny<=hist->GetNbinsY(); ++ibiny )
          for( Int_t ibinx=1; ibinx<=hist->GetNbinsX(); ++ibinx )
            {
             int globalbin = (ibiny-1)*hist->GetNbinsX() + ibinx - 1;
             m_HistoContents.set_fraction(globalbin,
  temp_HistoContents[globalbin]/integral);
  //          if(ihist>0)
  if(m_HistoContents.get_fraction(ihist-1)==m_HistoContents.get_fraction(ihist))
  {
  //            ATH_MSG_INFO("Skip fbin="<<ihist<<"
  fx="<<m_HistoBorders.GetBinLowEdge(ihist)<<"
  frac="<<m_HistoContents.get_fraction(ihist));
  //            --ihist;
  //            ATH_MSG_INFO("     fbin="<<ihist<<"
  fx="<<m_HistoBorders.GetBinLowEdge(ihist)<<"
  frac="<<m_HistoContents.get_fraction(ihist));
  //          }

  //          ATH_MSG_INFO("bin="<<ibin+1<<"
  fbin="<<ihist<<"/"<<m_HistoBorders.get_nbins()<<"
  x=["<<hist->GetXaxis()->GetBinLowEdge(ibin+1)<<","<<hist->GetXaxis()->GetBinUpEdge(ibin+1)<<"]
  fx="<<m_HistoBorders.GetBinLowEdge(ihist)<<"
  int="<<temp_HistoContents[ibin]/integral<<"
  frac="<<m_HistoContents.get_fraction(ihist));

  //          ++ihist;
                     }
  //        if(ihist!=nbins) {
  //          set_nbins(ihist);
  //
  m_HistoBorders.SetMinMax(hist->GetXaxis()->GetBinLowEdge(first+1),hist->GetXaxis()->GetBinUpEdge(last+1));
  //        }
          }//for ibinx
          }//for ibiny
  */
  // using TFCS2DFunction::rnd_to_fct;

  /// Function gets random number rnd in the range [0,1) as argument
  /// and returns function value according to a histogram distribution

  using TFCS2DFunction::rnd_to_fct;

  virtual void rnd_to_fct(float& valuex, float& valuey, float rnd) const
  {
    if (m_HistoContents.get_nbins() == 0) {
      valuex = 0.0;
      valuey = 0.0;
    }
    Trandom residual_rnd;
    size_t ibin = m_HistoContents.get_bin(rnd, residual_rnd);
    int ibiny = ibin % m_HistoBordersy.get_nbins();
    int ibinx = ibin / m_HistoBordersy.get_nbins();
    valuex = m_HistoBordersx.position(ibinx, 0.5);
    valuey = m_HistoBordersy.position(ibiny, residual_rnd);
    // valuex = 2.0+ 0.3*rnd;
    // valuey = 1200.0 + rnd*500.0;
  }

  virtual void rnd_to_fct(float& valuex,
                          float& valuey,
                          float rnd,
                          float rnd2) const
  {
    if (m_HistoContents.get_nbins() == 0) {
      valuex = 0.0;
      valuey = 0.0;
    }
    Trandom residual_rnd;
    size_t ibin = m_HistoContents.get_bin(rnd, residual_rnd);
    int ibiny = ibin
        % m_HistoBordersy
              .get_nbins();  // printf("%d ", m_HistoBordersy.get_nbins() );
    int ibinx = ibin
        / m_HistoBordersy
              .get_nbins();  // printf("%d ", m_HistoBordersx.get_nbins() );
                             // getchar();
    valuex = m_HistoBordersx.position(ibinx, rnd2);
    valuey = m_HistoBordersy.position(ibiny, residual_rnd);
    // valuex = 2.0+ rnd;
    // valuey = 1200.0 + rnd2*500.0;
  }
  /*
  virtual void rnd_to_fct(float& valuex,float& valuey,float rnd, float rnd2)
  const { if(m_HistoContents.get_nbins()==0) { valuex = 0.0; valuey = 0.0; }
        Trandom residual_rnd;
        size_t ibin=m_HistoContents.get_bin(rnd,residual_rnd);
        int ibiny = ibin / m_HistoBordersx.get_nbins();     //printf("%d ",
  m_HistoBordersy.get_nbins() ); int ibinx = ibin % m_HistoBordersx.get_nbins();
  //printf("%d ", m_HistoBordersx.get_nbins() ); getchar(); valuex =
  m_HistoBordersx.position(ibinx, rnd2); valuey =
  m_HistoBordersy.position(ibiny,residual_rnd);
          //valuex = 2.0+ rnd;
          //valuey = 1200.0 + rnd2*500.0;
  }*/

  inline const Txvec& get_HistoBordersx() const { return m_HistoBordersx; };
  inline const Tyvec& get_HistoBordersy() const { return m_HistoBordersy; };

  /////inline void get_XYZfromGlobal(int &ix, int &iy, int &iz, int global)
  /// const {hist;};

  inline Txvec& get_HistoBordersx() { return m_HistoBordersx; };
  inline Tyvec& get_HistoBordersy() { return m_HistoBordersy; };

  inline const TFCS1DFunction_HistogramContent<Tz, Trandom>& get_HistoContents()
      const
  {
    return m_HistoContents;
  };
  inline TFCS1DFunction_HistogramContent<Tz, Trandom>& get_HistoContents()
  {
    return m_HistoContents;
  };

protected:
  Txvec m_HistoBordersx;
  Tyvec m_HistoBordersy;
  TFCS1DFunction_HistogramContent<Tz, Trandom> m_HistoContents;

private:
  ClassDef(TFCS2DFunctionTemplateHistogram,
           1)  // TFCS2DFunctionTemplateHistogram
};

class TFCS2DFunctionInt8Int8Int8Histogram
    : public TFCS2DFunctionTemplateHistogram<
          TFCS1DFunction_HistogramInt8BinEdges,
          TFCS1DFunction_HistogramInt8BinEdges,
          uint8_t,
          float>
{
public:
  TFCS2DFunctionInt8Int8Int8Histogram(TH2* h = nullptr)
      : TFCS2DFunctionTemplateHistogram<TFCS1DFunction_HistogramInt8BinEdges,
                                        TFCS1DFunction_HistogramInt8BinEdges,
                                        uint8_t,
                                        float>(h) {};

  ClassDef(TFCS2DFunctionInt8Int8Int8Histogram,
           1)  // TFCS2DFunctionInt8Int8Int8Histogram
};

class TFCS2DFunctionInt8Int8Int16Histogram
    : public TFCS2DFunctionTemplateHistogram<
          TFCS1DFunction_HistogramInt8BinEdges,
          TFCS1DFunction_HistogramInt8BinEdges,
          uint16_t,
          float>
{
public:
  TFCS2DFunctionInt8Int8Int16Histogram(TH2* h = nullptr)
      : TFCS2DFunctionTemplateHistogram<TFCS1DFunction_HistogramInt8BinEdges,
                                        TFCS1DFunction_HistogramInt8BinEdges,
                                        uint16_t,
                                        float>(h) {};

  ClassDef(TFCS2DFunctionInt8Int8Int16Histogram,
           1)  // TFCS2DFunctionInt8Int8Int16Histogram
};

class TFCS2DFunctionInt8Int8Int32Histogram
    : public TFCS2DFunctionTemplateHistogram<
          TFCS1DFunction_HistogramInt8BinEdges,
          TFCS1DFunction_HistogramInt8BinEdges,
          uint32_t,
          float>
{
public:
  TFCS2DFunctionInt8Int8Int32Histogram(TH2* h = nullptr)
      : TFCS2DFunctionTemplateHistogram<TFCS1DFunction_HistogramInt8BinEdges,
                                        TFCS1DFunction_HistogramInt8BinEdges,
                                        uint32_t,
                                        float>(h) {};

  ClassDef(TFCS2DFunctionInt8Int8Int32Histogram,
           1)  // TFCS2DFunctionInt8Int8Int32Histogram
};
class TFCS2DFunctionInt8Int16Int8Histogram
    : public TFCS2DFunctionTemplateHistogram<
          TFCS1DFunction_HistogramInt8BinEdges,
          TFCS1DFunction_HistogramInt16BinEdges,
          uint8_t,
          float>
{
public:
  TFCS2DFunctionInt8Int16Int8Histogram(TH2* h = nullptr)
      : TFCS2DFunctionTemplateHistogram<TFCS1DFunction_HistogramInt8BinEdges,
                                        TFCS1DFunction_HistogramInt16BinEdges,
                                        uint8_t,
                                        float>(h) {};

  ClassDef(TFCS2DFunctionInt8Int16Int8Histogram,
           1)  // TFCS2DFunctionInt8Int16Int8Histogram
};

class TFCS2DFunctionInt8Int16Int16Histogram
    : public TFCS2DFunctionTemplateHistogram<
          TFCS1DFunction_HistogramInt8BinEdges,
          TFCS1DFunction_HistogramInt16BinEdges,
          uint16_t,
          float>
{
public:
  TFCS2DFunctionInt8Int16Int16Histogram(TH2* h = nullptr)
      : TFCS2DFunctionTemplateHistogram<TFCS1DFunction_HistogramInt8BinEdges,
                                        TFCS1DFunction_HistogramInt16BinEdges,
                                        uint16_t,
                                        float>(h) {};

  ClassDef(TFCS2DFunctionInt8Int16Int16Histogram,
           1)  // TFCS2DFunctionInt8Int16Int16Histogram
};

class TFCS2DFunctionInt8Int16Int32Histogram
    : public TFCS2DFunctionTemplateHistogram<
          TFCS1DFunction_HistogramInt8BinEdges,
          TFCS1DFunction_HistogramInt16BinEdges,
          uint32_t,
          float>
{
public:
  TFCS2DFunctionInt8Int16Int32Histogram(TH2* h = nullptr)
      : TFCS2DFunctionTemplateHistogram<TFCS1DFunction_HistogramInt8BinEdges,
                                        TFCS1DFunction_HistogramInt16BinEdges,
                                        uint32_t,
                                        float>(h) {};

  ClassDef(TFCS2DFunctionInt8Int16Int32Histogram,
           1)  // TFCS2DFunctionInt8Int16Int32Histogram
};
class TFCS2DFunctionInt8Int32Int8Histogram
    : public TFCS2DFunctionTemplateHistogram<
          TFCS1DFunction_HistogramInt8BinEdges,
          TFCS1DFunction_HistogramInt32BinEdges,
          uint8_t,
          float>
{
public:
  TFCS2DFunctionInt8Int32Int8Histogram(TH2* h = nullptr)
      : TFCS2DFunctionTemplateHistogram<TFCS1DFunction_HistogramInt8BinEdges,
                                        TFCS1DFunction_HistogramInt32BinEdges,
                                        uint8_t,
                                        float>(h) {};

  ClassDef(TFCS2DFunctionInt8Int32Int8Histogram,
           1)  // TFCS2DFunctionInt8Int32Int8Histogram
};
class TFCS2DFunctionInt8Int32Int16Histogram
    : public TFCS2DFunctionTemplateHistogram<
          TFCS1DFunction_HistogramInt8BinEdges,
          TFCS1DFunction_HistogramInt32BinEdges,
          uint16_t,
          float>
{
public:
  TFCS2DFunctionInt8Int32Int16Histogram(TH2* h = nullptr)
      : TFCS2DFunctionTemplateHistogram<TFCS1DFunction_HistogramInt8BinEdges,
                                        TFCS1DFunction_HistogramInt32BinEdges,
                                        uint16_t,
                                        float>(h) {};

  ClassDef(TFCS2DFunctionInt8Int32Int16Histogram,
           1)  // TFCS2DFunctionInt8Int32Int16Histogram
};
class TFCS2DFunctionInt8Int32Int32Histogram
    : public TFCS2DFunctionTemplateHistogram<
          TFCS1DFunction_HistogramInt8BinEdges,
          TFCS1DFunction_HistogramInt32BinEdges,
          uint32_t,
          float>
{
public:
  TFCS2DFunctionInt8Int32Int32Histogram(TH2* h = nullptr)
      : TFCS2DFunctionTemplateHistogram<TFCS1DFunction_HistogramInt8BinEdges,
                                        TFCS1DFunction_HistogramInt32BinEdges,
                                        uint32_t,
                                        float>(h) {};

  ClassDef(TFCS2DFunctionInt8Int32Int32Histogram,
           1)  // TFCS2DFunctionInt8Int32Int32Histogram
};

class TFCS2DFunctionInt16Int8Int8Histogram
    : public TFCS2DFunctionTemplateHistogram<
          TFCS1DFunction_HistogramInt16BinEdges,
          TFCS1DFunction_HistogramInt8BinEdges,
          uint8_t,
          float>
{
public:
  TFCS2DFunctionInt16Int8Int8Histogram(TH2* h = nullptr)
      : TFCS2DFunctionTemplateHistogram<TFCS1DFunction_HistogramInt16BinEdges,
                                        TFCS1DFunction_HistogramInt8BinEdges,
                                        uint8_t,
                                        float>(h) {};

  ClassDef(TFCS2DFunctionInt16Int8Int8Histogram,
           1)  // TFCS2DFunctionInt16Int8Int8Histogram
};

class TFCS2DFunctionInt16Int8Int16Histogram
    : public TFCS2DFunctionTemplateHistogram<
          TFCS1DFunction_HistogramInt16BinEdges,
          TFCS1DFunction_HistogramInt8BinEdges,
          uint16_t,
          float>
{
public:
  TFCS2DFunctionInt16Int8Int16Histogram(TH2* h = nullptr)
      : TFCS2DFunctionTemplateHistogram<TFCS1DFunction_HistogramInt16BinEdges,
                                        TFCS1DFunction_HistogramInt8BinEdges,
                                        uint16_t,
                                        float>(h) {};

  ClassDef(TFCS2DFunctionInt16Int8Int16Histogram,
           1)  // TFCS2DFunctionInt16Int8Int16Histogram
};

class TFCS2DFunctionInt16Int8Int32Histogram
    : public TFCS2DFunctionTemplateHistogram<
          TFCS1DFunction_HistogramInt16BinEdges,
          TFCS1DFunction_HistogramInt8BinEdges,
          uint32_t,
          float>
{
public:
  TFCS2DFunctionInt16Int8Int32Histogram(TH2* h = nullptr)
      : TFCS2DFunctionTemplateHistogram<TFCS1DFunction_HistogramInt16BinEdges,
                                        TFCS1DFunction_HistogramInt8BinEdges,
                                        uint32_t,
                                        float>(h) {};

  ClassDef(TFCS2DFunctionInt16Int8Int32Histogram,
           1)  // TFCS2DFunctionInt16Int8Int32Histogram
};
class TFCS2DFunctionInt16Int16Int8Histogram
    : public TFCS2DFunctionTemplateHistogram<
          TFCS1DFunction_HistogramInt16BinEdges,
          TFCS1DFunction_HistogramInt16BinEdges,
          uint8_t,
          float>
{
public:
  TFCS2DFunctionInt16Int16Int8Histogram(TH2* h = nullptr)
      : TFCS2DFunctionTemplateHistogram<TFCS1DFunction_HistogramInt16BinEdges,
                                        TFCS1DFunction_HistogramInt16BinEdges,
                                        uint8_t,
                                        float>(h) {};

  ClassDef(TFCS2DFunctionInt16Int16Int8Histogram,
           1)  // TFCS2DFunctionInt16Int16Int8Histogram
};

class TFCS2DFunctionInt16Int16Int16Histogram
    : public TFCS2DFunctionTemplateHistogram<
          TFCS1DFunction_HistogramInt16BinEdges,
          TFCS1DFunction_HistogramInt16BinEdges,
          uint16_t,
          float>
{
public:
  TFCS2DFunctionInt16Int16Int16Histogram(TH2* h = nullptr)
      : TFCS2DFunctionTemplateHistogram<TFCS1DFunction_HistogramInt16BinEdges,
                                        TFCS1DFunction_HistogramInt16BinEdges,
                                        uint16_t,
                                        float>(h) {};

  ClassDef(TFCS2DFunctionInt16Int16Int16Histogram,
           1)  // TFCS2DFunctionInt16Int16Int16Histogram
};

class TFCS2DFunctionInt16Int16Int32Histogram
    : public TFCS2DFunctionTemplateHistogram<
          TFCS1DFunction_HistogramInt16BinEdges,
          TFCS1DFunction_HistogramInt16BinEdges,
          uint32_t,
          float>
{
public:
  TFCS2DFunctionInt16Int16Int32Histogram(TH2* h = nullptr)
      : TFCS2DFunctionTemplateHistogram<TFCS1DFunction_HistogramInt16BinEdges,
                                        TFCS1DFunction_HistogramInt16BinEdges,
                                        uint32_t,
                                        float>(h) {};

  ClassDef(TFCS2DFunctionInt16Int16Int32Histogram,
           1)  // TFCS2DFunctionInt16Int16Int32Histogram
};
class TFCS2DFunctionInt16Int32Int8Histogram
    : public TFCS2DFunctionTemplateHistogram<
          TFCS1DFunction_HistogramInt16BinEdges,
          TFCS1DFunction_HistogramInt32BinEdges,
          uint8_t,
          float>
{
public:
  TFCS2DFunctionInt16Int32Int8Histogram(TH2* h = nullptr)
      : TFCS2DFunctionTemplateHistogram<TFCS1DFunction_HistogramInt16BinEdges,
                                        TFCS1DFunction_HistogramInt32BinEdges,
                                        uint8_t,
                                        float>(h) {};

  ClassDef(TFCS2DFunctionInt16Int32Int8Histogram,
           1)  // TFCS2DFunctionInt16Int32Int8Histogram
};
class TFCS2DFunctionInt16Int32Int16Histogram
    : public TFCS2DFunctionTemplateHistogram<
          TFCS1DFunction_HistogramInt16BinEdges,
          TFCS1DFunction_HistogramInt32BinEdges,
          uint16_t,
          float>
{
public:
  TFCS2DFunctionInt16Int32Int16Histogram(TH2* h = nullptr)
      : TFCS2DFunctionTemplateHistogram<TFCS1DFunction_HistogramInt16BinEdges,
                                        TFCS1DFunction_HistogramInt32BinEdges,
                                        uint16_t,
                                        float>(h) {};

  ClassDef(TFCS2DFunctionInt16Int32Int16Histogram,
           1)  // TFCS2DFunctionInt16Int32Int16Histogram
};
class TFCS2DFunctionInt16Int32Int32Histogram
    : public TFCS2DFunctionTemplateHistogram<
          TFCS1DFunction_HistogramInt16BinEdges,
          TFCS1DFunction_HistogramInt32BinEdges,
          uint32_t,
          float>
{
public:
  TFCS2DFunctionInt16Int32Int32Histogram(TH2* h = nullptr)
      : TFCS2DFunctionTemplateHistogram<TFCS1DFunction_HistogramInt16BinEdges,
                                        TFCS1DFunction_HistogramInt32BinEdges,
                                        uint32_t,
                                        float>(h) {};

  ClassDef(TFCS2DFunctionInt16Int32Int32Histogram,
           1)  // TFCS2DFunctionInt16Int32Int32Histogram
};

class TFCS2DFunctionInt32Int8Int8Histogram
    : public TFCS2DFunctionTemplateHistogram<
          TFCS1DFunction_HistogramInt32BinEdges,
          TFCS1DFunction_HistogramInt8BinEdges,
          uint8_t,
          float>
{
public:
  TFCS2DFunctionInt32Int8Int8Histogram(TH2* h = nullptr)
      : TFCS2DFunctionTemplateHistogram<TFCS1DFunction_HistogramInt32BinEdges,
                                        TFCS1DFunction_HistogramInt8BinEdges,
                                        uint8_t,
                                        float>(h) {};

  ClassDef(TFCS2DFunctionInt32Int8Int8Histogram,
           1)  // TFCS2DFunctionInt32Int8Int8Histogram
};

class TFCS2DFunctionInt32Int8Int16Histogram
    : public TFCS2DFunctionTemplateHistogram<
          TFCS1DFunction_HistogramInt32BinEdges,
          TFCS1DFunction_HistogramInt8BinEdges,
          uint16_t,
          float>
{
public:
  TFCS2DFunctionInt32Int8Int16Histogram(TH2* h = nullptr)
      : TFCS2DFunctionTemplateHistogram<TFCS1DFunction_HistogramInt32BinEdges,
                                        TFCS1DFunction_HistogramInt8BinEdges,
                                        uint16_t,
                                        float>(h) {};

  ClassDef(TFCS2DFunctionInt32Int8Int16Histogram,
           1)  // TFCS2DFunctionInt32Int8Int16Histogram
};

class TFCS2DFunctionInt32Int8Int32Histogram
    : public TFCS2DFunctionTemplateHistogram<
          TFCS1DFunction_HistogramInt32BinEdges,
          TFCS1DFunction_HistogramInt8BinEdges,
          uint32_t,
          float>
{
public:
  TFCS2DFunctionInt32Int8Int32Histogram(TH2* h = nullptr)
      : TFCS2DFunctionTemplateHistogram<TFCS1DFunction_HistogramInt32BinEdges,
                                        TFCS1DFunction_HistogramInt8BinEdges,
                                        uint32_t,
                                        float>(h) {};

  ClassDef(TFCS2DFunctionInt32Int8Int32Histogram,
           1)  // TFCS2DFunctionInt32Int8Int32Histogram
};
class TFCS2DFunctionInt32Int16Int8Histogram
    : public TFCS2DFunctionTemplateHistogram<
          TFCS1DFunction_HistogramInt32BinEdges,
          TFCS1DFunction_HistogramInt16BinEdges,
          uint8_t,
          float>
{
public:
  TFCS2DFunctionInt32Int16Int8Histogram(TH2* h = nullptr)
      : TFCS2DFunctionTemplateHistogram<TFCS1DFunction_HistogramInt32BinEdges,
                                        TFCS1DFunction_HistogramInt16BinEdges,
                                        uint8_t,
                                        float>(h) {};

  ClassDef(TFCS2DFunctionInt32Int16Int8Histogram,
           1)  // TFCS2DFunctionInt32Int16Int8Histogram
};

class TFCS2DFunctionInt32Int16Int16Histogram
    : public TFCS2DFunctionTemplateHistogram<
          TFCS1DFunction_HistogramInt32BinEdges,
          TFCS1DFunction_HistogramInt16BinEdges,
          uint16_t,
          float>
{
public:
  TFCS2DFunctionInt32Int16Int16Histogram(TH2* h = nullptr)
      : TFCS2DFunctionTemplateHistogram<TFCS1DFunction_HistogramInt32BinEdges,
                                        TFCS1DFunction_HistogramInt16BinEdges,
                                        uint16_t,
                                        float>(h) {};

  ClassDef(TFCS2DFunctionInt32Int16Int16Histogram,
           1)  // TFCS2DFunctionInt32Int16Int16Histogram
};

class TFCS2DFunctionInt32Int16Int32Histogram
    : public TFCS2DFunctionTemplateHistogram<
          TFCS1DFunction_HistogramInt32BinEdges,
          TFCS1DFunction_HistogramInt16BinEdges,
          uint32_t,
          float>
{
public:
  TFCS2DFunctionInt32Int16Int32Histogram(TH2* h = nullptr)
      : TFCS2DFunctionTemplateHistogram<TFCS1DFunction_HistogramInt32BinEdges,
                                        TFCS1DFunction_HistogramInt16BinEdges,
                                        uint32_t,
                                        float>(h) {};

  ClassDef(TFCS2DFunctionInt32Int16Int32Histogram,
           1)  // TFCS2DFunctionInt32Int16Int32Histogram
};
class TFCS2DFunctionInt32Int32Int8Histogram
    : public TFCS2DFunctionTemplateHistogram<
          TFCS1DFunction_HistogramInt32BinEdges,
          TFCS1DFunction_HistogramInt32BinEdges,
          uint8_t,
          float>
{
public:
  TFCS2DFunctionInt32Int32Int8Histogram(TH2* h = nullptr)
      : TFCS2DFunctionTemplateHistogram<TFCS1DFunction_HistogramInt32BinEdges,
                                        TFCS1DFunction_HistogramInt32BinEdges,
                                        uint8_t,
                                        float>(h) {};

  ClassDef(TFCS2DFunctionInt32Int32Int8Histogram,
           1)  // TFCS2DFunctionInt32Int32Int8Histogram
};
class TFCS2DFunctionInt32Int32Int16Histogram
    : public TFCS2DFunctionTemplateHistogram<
          TFCS1DFunction_HistogramInt32BinEdges,
          TFCS1DFunction_HistogramInt32BinEdges,
          uint16_t,
          float>
{
public:
  TFCS2DFunctionInt32Int32Int16Histogram(TH2* h = nullptr)
      : TFCS2DFunctionTemplateHistogram<TFCS1DFunction_HistogramInt32BinEdges,
                                        TFCS1DFunction_HistogramInt32BinEdges,
                                        uint16_t,
                                        float>(h) {};

  ClassDef(TFCS2DFunctionInt32Int32Int16Histogram,
           1)  // TFCS2DFunctionInt32Int32Int16Histogram
};
class TFCS2DFunctionInt32Int32Int32Histogram
    : public TFCS2DFunctionTemplateHistogram<
          TFCS1DFunction_HistogramInt32BinEdges,
          TFCS1DFunction_HistogramInt32BinEdges,
          uint32_t,
          float>
{
public:
  TFCS2DFunctionInt32Int32Int32Histogram(TH2* h = nullptr)
      : TFCS2DFunctionTemplateHistogram<TFCS1DFunction_HistogramInt32BinEdges,
                                        TFCS1DFunction_HistogramInt32BinEdges,
                                        uint32_t,
                                        float>(h) {};

  ClassDef(TFCS2DFunctionInt32Int32Int32Histogram,
           1)  // TFCS2DFunctionInt32Int32Int32Histogram
};

// clang-format off
#if defined(__ROOTCLING__) && defined(__FastCaloSimStandAlone__)
#  pragma link C++ class TFCS2DFunctionTemplateHistogram \
          < TFCS1DFunction_HistogramInt8BinEdges, \
      TFCS1DFunction_HistogramInt8BinEdges, uint8_t, float> + ;
#  pragma link C++ class TFCS2DFunctionTemplateHistogram \
          < TFCS1DFunction_HistogramInt8BinEdges, \
      TFCS1DFunction_HistogramInt8BinEdges, uint16_t, float> + ;
#  pragma link C++ class TFCS2DFunctionTemplateHistogram \
          < TFCS1DFunction_HistogramInt8BinEdges, \
      TFCS1DFunction_HistogramInt8BinEdges, uint32_t, float> + ;
#  pragma link C++ class TFCS2DFunctionTemplateHistogram \
          < TFCS1DFunction_HistogramInt8BinEdges, \
      TFCS1DFunction_HistogramInt16BinEdges, uint8_t, float> + ;
#  pragma link C++ class TFCS2DFunctionTemplateHistogram \
          < TFCS1DFunction_HistogramInt8BinEdges, \
      TFCS1DFunction_HistogramInt16BinEdges, uint16_t, float> + ;
#  pragma link C++ class TFCS2DFunctionTemplateHistogram \
          < TFCS1DFunction_HistogramInt8BinEdges, \
      TFCS1DFunction_HistogramInt16BinEdges, uint32_t, float> + ;
#  pragma link C++ class TFCS2DFunctionTemplateHistogram \
          < TFCS1DFunction_HistogramInt8BinEdges, \
      TFCS1DFunction_HistogramInt32BinEdges, uint8_t, float> + ;
#  pragma link C++ class TFCS2DFunctionTemplateHistogram \
          < TFCS1DFunction_HistogramInt8BinEdges, \
      TFCS1DFunction_HistogramInt32BinEdges, uint16_t, float> + ;
#  pragma link C++ class TFCS2DFunctionTemplateHistogram \
          < TFCS1DFunction_HistogramInt8BinEdges, \
      TFCS1DFunction_HistogramInt32BinEdges, uint32_t, float> + ;
#  pragma link C++ class TFCS2DFunctionTemplateHistogram \
          < TFCS1DFunction_HistogramInt16BinEdges, \
      TFCS1DFunction_HistogramInt8BinEdges, uint8_t, float> + ;
#  pragma link C++ class TFCS2DFunctionTemplateHistogram \
          < TFCS1DFunction_HistogramInt16BinEdges, \
      TFCS1DFunction_HistogramInt8BinEdges, uint16_t, float> + ;
#  pragma link C++ class TFCS2DFunctionTemplateHistogram \
          < TFCS1DFunction_HistogramInt16BinEdges, \
      TFCS1DFunction_HistogramInt8BinEdges, uint32_t, float> + ;
#  pragma link C++ class TFCS2DFunctionTemplateHistogram \
          < TFCS1DFunction_HistogramInt16BinEdges, \
      TFCS1DFunction_HistogramInt16BinEdges, uint8_t, float> + ;
#  pragma link C++ class TFCS2DFunctionTemplateHistogram \
          < TFCS1DFunction_HistogramInt16BinEdges, \
      TFCS1DFunction_HistogramInt16BinEdges, uint16_t, float> + ;
#  pragma link C++ class TFCS2DFunctionTemplateHistogram \
          < TFCS1DFunction_HistogramInt16BinEdges, \
      TFCS1DFunction_HistogramInt16BinEdges, uint32_t, float> + ;
#  pragma link C++ class TFCS2DFunctionTemplateHistogram \
          < TFCS1DFunction_HistogramInt16BinEdges, \
      TFCS1DFunction_HistogramInt32BinEdges, uint8_t, float> + ;
#  pragma link C++ class TFCS2DFunctionTemplateHistogram \
          < TFCS1DFunction_HistogramInt16BinEdges, \
      TFCS1DFunction_HistogramInt32BinEdges, uint16_t, float> + ;
#  pragma link C++ class TFCS2DFunctionTemplateHistogram \
          < TFCS1DFunction_HistogramInt16BinEdges, \
      TFCS1DFunction_HistogramInt32BinEdges, uint32_t, float> + ;
#  pragma link C++ class TFCS2DFunctionTemplateHistogram \
          < TFCS1DFunction_HistogramInt32BinEdges, \
      TFCS1DFunction_HistogramInt8BinEdges, uint8_t, float> + ;
#  pragma link C++ class TFCS2DFunctionTemplateHistogram \
          < TFCS1DFunction_HistogramInt32BinEdges, \
      TFCS1DFunction_HistogramInt8BinEdges, uint16_t, float> + ;
#  pragma link C++ class TFCS2DFunctionTemplateHistogram \
          < TFCS1DFunction_HistogramInt32BinEdges, \
      TFCS1DFunction_HistogramInt8BinEdges, uint32_t, float> + ;
#  pragma link C++ class TFCS2DFunctionTemplateHistogram \
          < TFCS1DFunction_HistogramInt32BinEdges, \
      TFCS1DFunction_HistogramInt16BinEdges, uint8_t, float> + ;
#  pragma link C++ class TFCS2DFunctionTemplateHistogram \
          < TFCS1DFunction_HistogramInt32BinEdges, \
      TFCS1DFunction_HistogramInt16BinEdges, uint16_t, float> + ;
#  pragma link C++ class TFCS2DFunctionTemplateHistogram \
          < TFCS1DFunction_HistogramInt32BinEdges, \
      TFCS1DFunction_HistogramInt16BinEdges, uint32_t, float> + ;
#  pragma link C++ class TFCS2DFunctionTemplateHistogram \
          < TFCS1DFunction_HistogramInt32BinEdges, \
      TFCS1DFunction_HistogramInt32BinEdges, uint8_t, float> + ;
#  pragma link C++ class TFCS2DFunctionTemplateHistogram \
          < TFCS1DFunction_HistogramInt32BinEdges, \
      TFCS1DFunction_HistogramInt32BinEdges, uint16_t, float> + ;
#  pragma link C++ class TFCS2DFunctionTemplateHistogram \
          < TFCS1DFunction_HistogramInt32BinEdges, \
      TFCS1DFunction_HistogramInt32BinEdges, uint32_t, float> + ;

#  pragma link C++ class TFCS2DFunctionInt8Int8Int8Histogram + ;
#  pragma link C++ class TFCS2DFunctionInt8Int8Int16Histogram + ;
#  pragma link C++ class TFCS2DFunctionInt8Int8Int32Histogram + ;
#  pragma link C++ class TFCS2DFunctionInt8Int16Int8Histogram + ;
#  pragma link C++ class TFCS2DFunctionInt8Int16Int16Histogram + ;
#  pragma link C++ class TFCS2DFunctionInt8Int16Int32Histogram + ;
#  pragma link C++ class TFCS2DFunctionInt8Int32Int8Histogram + ;
#  pragma link C++ class TFCS2DFunctionInt8Int32Int16Histogram + ;
#  pragma link C++ class TFCS2DFunctionInt8Int32Int32Histogram + ;

#  pragma link C++ class TFCS2DFunctionInt16Int8Int8Histogram + ;
#  pragma link C++ class TFCS2DFunctionInt16Int8Int16Histogram + ;
#  pragma link C++ class TFCS2DFunctionInt16Int8Int32Histogram + ;
#  pragma link C++ class TFCS2DFunctionInt16Int16Int8Histogram + ;
#  pragma link C++ class TFCS2DFunctionInt16Int16Int16Histogram + ;
#  pragma link C++ class TFCS2DFunctionInt16Int16Int32Histogram + ;
#  pragma link C++ class TFCS2DFunctionInt16Int32Int8Histogram + ;
#  pragma link C++ class TFCS2DFunctionInt16Int32Int16Histogram + ;
#  pragma link C++ class TFCS2DFunctionInt16Int32Int32Histogram + ;

#  pragma link C++ class TFCS2DFunctionInt32Int8Int8Histogram + ;
#  pragma link C++ class TFCS2DFunctionInt32Int8Int16Histogram + ;
#  pragma link C++ class TFCS2DFunctionInt32Int8Int32Histogram + ;
#  pragma link C++ class TFCS2DFunctionInt32Int16Int8Histogram + ;
#  pragma link C++ class TFCS2DFunctionInt32Int16Int16Histogram + ;
#  pragma link C++ class TFCS2DFunctionInt32Int16Int32Histogram + ;
#  pragma link C++ class TFCS2DFunctionInt32Int32Int8Histogram + ;
#  pragma link C++ class TFCS2DFunctionInt32Int32Int16Histogram + ;
#  pragma link C++ class TFCS2DFunctionInt32Int32Int32Histogram + ;

#endif
// clang-format on

#endif
