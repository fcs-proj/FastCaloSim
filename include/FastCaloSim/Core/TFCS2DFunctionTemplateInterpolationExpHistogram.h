// Copyright (c) 2024 CERN for the benefit of the FastCaloSim project

#ifndef TFCS2DFunctionTemplateInterpolationExpHistogram_h
#define TFCS2DFunctionTemplateInterpolationExpHistogram_h

#include <iostream>

#include "FastCaloSim/Core/TFCS2DFunctionTemplateHistogram.h"
#include "TBuffer.h"
#include "TFCS2DFunction.h"
#include "TH1.h"
#include "TH2.h"

namespace FastCaloSim::Core
{
template<typename Txvec, typename Tyvec, typename Tz, typename Trandom = float>
class TFCS2DFunctionTemplateInterpolationExpHistogram
    : public TFCS2DFunctionTemplateHistogram<Txvec, Tyvec, Tz, Trandom>
{
public:
  TFCS2DFunctionTemplateInterpolationExpHistogram(TH2* hist = nullptr)
      : TFCS2DFunctionTemplateHistogram<Txvec, Tyvec, Tz, Trandom>(hist) {};

  using TFCS2DFunctionTemplateHistogram<Txvec, Tyvec, Tz, Trandom>::rnd_to_fct;
  using TFCS2DFunctionTemplateHistogram<Txvec, Tyvec, Tz, Trandom>::get_nbins;
  using TFCS2DFunctionTemplateHistogram<Txvec, Tyvec, Tz, Trandom>::
      m_HistoContents;
  using TFCS2DFunctionTemplateHistogram<Txvec, Tyvec, Tz, Trandom>::
      m_HistoBordersx;
  using TFCS2DFunctionTemplateHistogram<Txvec, Tyvec, Tz, Trandom>::
      m_HistoBordersy;

  /// Function gets random number rnd in the range [0,1) as argument
  /// and returns function value according to a histogram distribution.
  /// An exponential interpolation is done within each bin

  virtual void rnd_to_fct(float& valuex,
                          float& valuey,
                          float rnd,
                          float rnd2) const
  {
    size_t nbins = get_nbins();
    if (nbins == 0) {
      valuex = 0.0;
      valuey = 0.0;
      return;
    }
    Trandom residual_rnd;
    size_t ibin = m_HistoContents.get_bin(rnd, residual_rnd);

    Trandom frac = m_HistoContents.get_fraction(ibin);

    size_t ibiny = ibin
        % m_HistoBordersy
              .get_nbins();  // printf("%d ", m_HistoBordersy.get_nbins() );
    size_t ibinx = ibin / m_HistoBordersy.get_nbins();
    // size_t nbinsx = m_HistoBordersx.get_nbins();
    size_t nbinsy = m_HistoBordersy.get_nbins();

    Trandom dfrac = 0.0;
    auto l = m_HistoBordersy.GetBinLength(ibiny);
    Trandom dfracprev = 0.0;
    Trandom ldfrac = 0.0;
    auto lprev = l;
    if (ibiny > 0) {
      Trandom fracprev = m_HistoContents.get_fraction(ibin - 1);
      dfrac = frac - fracprev;
      if (ibiny > 1) {
        dfracprev = fracprev - m_HistoContents.get_fraction(ibin - 2);
      } else {
        dfracprev = fracprev;
        if (ibinx > 0)
          dfracprev -= m_HistoContents.get_fraction(ibin - 2);
      }
      lprev = m_HistoBordersy.GetBinLength(ibiny - 1);
    } else {
      dfrac = frac;
      if (ibinx > 0)
        dfrac -= m_HistoContents.get_fraction(ibin - 1);
      dfracprev = dfrac;
    }

    if (dfrac > 0.0 && dfracprev > 0.0) {
      ldfrac = log(dfrac / l) - log(dfracprev / lprev);
    }

    Trandom dfracnext = 0.0;
    Trandom ldfracnext = 0.0;
    auto lnext = l;
    if (ibiny < nbinsy - 1) {
      dfracnext = m_HistoContents.get_fraction(ibin + 1) - frac;
      lnext = m_HistoBordersy.GetBinLength(ibiny + 1);
    } else {
      dfracnext = dfrac;
    }

    if (dfrac > 0.0 && dfracnext > 0.0) {
      ldfracnext = log(dfracnext / lnext) - log(dfrac / l);
    }

    Trandom beta = 0.0, betaprev = 0.0, betanext = 0.0;
    betanext = ldfracnext / (lnext + l) / 0.5;
    betaprev = ldfrac / (lprev + l) / 0.5;
    beta = 0.5 * (betaprev + betanext);
    // printf("betaprev=%13.5e betanext=%13.5e beta=%13.5e\n", betaprev,
    // betanext, beta); printf("l=%13.5e lnext=%13.5e\n", l, lnext);
    // printf("log=%13.5e lognext=%13.5e logprev=%13.5e\n", log(dfrac/l),
    // log(dfracnext/lnext), log(dfracprev/lprev) ); printf("ldfrac=%13.5e
    // ldfracnext=%13.5e\n", ldfrac, ldfracnext );

    valuey = m_HistoBordersy.position_exp(ibiny, beta, residual_rnd);
    valuex = m_HistoBordersx.position(ibinx, rnd2);
  }

  ClassDef(TFCS2DFunctionTemplateInterpolationExpHistogram,
           1)  // TFCS1DFunctionTemplateInterpolationExpHistogram
};

class TFCS2DFunctionInt8Int8Int8InterpolationExpHistogram
    : public TFCS2DFunctionTemplateInterpolationExpHistogram<
          TFCS1DFunction_HistogramInt8BinEdges,
          TFCS1DFunction_HistogramInt8BinEdges,
          uint8_t,
          float>
{
public:
  TFCS2DFunctionInt8Int8Int8InterpolationExpHistogram(TH2* h = nullptr)
      : TFCS2DFunctionTemplateInterpolationExpHistogram<
            TFCS1DFunction_HistogramInt8BinEdges,
            TFCS1DFunction_HistogramInt8BinEdges,
            uint8_t,
            float>(h) {};

  ClassDef(TFCS2DFunctionInt8Int8Int8InterpolationExpHistogram,
           1)  // TFCS2DFunctionInt8Int8Int8InterpolationExpHistogram
};
class TFCS2DFunctionInt8Int8Int16InterpolationExpHistogram
    : public TFCS2DFunctionTemplateInterpolationExpHistogram<
          TFCS1DFunction_HistogramInt8BinEdges,
          TFCS1DFunction_HistogramInt8BinEdges,
          uint16_t,
          float>
{
public:
  TFCS2DFunctionInt8Int8Int16InterpolationExpHistogram(TH2* h = nullptr)
      : TFCS2DFunctionTemplateInterpolationExpHistogram<
            TFCS1DFunction_HistogramInt8BinEdges,
            TFCS1DFunction_HistogramInt8BinEdges,
            uint16_t,
            float>(h) {};

  ClassDef(TFCS2DFunctionInt8Int8Int16InterpolationExpHistogram,
           1)  // TFCS2DFunctionInt8Int8Int16InterpolationExpHistogram
};
class TFCS2DFunctionInt8Int8Int32InterpolationExpHistogram
    : public TFCS2DFunctionTemplateInterpolationExpHistogram<
          TFCS1DFunction_HistogramInt8BinEdges,
          TFCS1DFunction_HistogramInt8BinEdges,
          uint32_t,
          float>
{
public:
  TFCS2DFunctionInt8Int8Int32InterpolationExpHistogram(TH2* h = nullptr)
      : TFCS2DFunctionTemplateInterpolationExpHistogram<
            TFCS1DFunction_HistogramInt8BinEdges,
            TFCS1DFunction_HistogramInt8BinEdges,
            uint32_t,
            float>(h) {};

  ClassDef(TFCS2DFunctionInt8Int8Int32InterpolationExpHistogram,
           1)  // TFCS2DFunctionInt8Int8Int32InterpolationExpHistogram
};

////

class TFCS2DFunctionInt8Int16Int8InterpolationExpHistogram
    : public TFCS2DFunctionTemplateInterpolationExpHistogram<
          TFCS1DFunction_HistogramInt8BinEdges,
          TFCS1DFunction_HistogramInt16BinEdges,
          uint8_t,
          float>
{
public:
  TFCS2DFunctionInt8Int16Int8InterpolationExpHistogram(TH2* h = nullptr)
      : TFCS2DFunctionTemplateInterpolationExpHistogram<
            TFCS1DFunction_HistogramInt8BinEdges,
            TFCS1DFunction_HistogramInt16BinEdges,
            uint8_t,
            float>(h) {};

  ClassDef(TFCS2DFunctionInt8Int16Int8InterpolationExpHistogram,
           1)  // TFCS2DFunctionInt8Int16Int8InterpolationExpHistogram
};
class TFCS2DFunctionInt8Int16Int16InterpolationExpHistogram
    : public TFCS2DFunctionTemplateInterpolationExpHistogram<
          TFCS1DFunction_HistogramInt8BinEdges,
          TFCS1DFunction_HistogramInt16BinEdges,
          uint16_t,
          float>
{
public:
  TFCS2DFunctionInt8Int16Int16InterpolationExpHistogram(TH2* h = nullptr)
      : TFCS2DFunctionTemplateInterpolationExpHistogram<
            TFCS1DFunction_HistogramInt8BinEdges,
            TFCS1DFunction_HistogramInt16BinEdges,
            uint16_t,
            float>(h) {};

  ClassDef(TFCS2DFunctionInt8Int16Int16InterpolationExpHistogram,
           1)  // TFCS2DFunctionInt8Int16Int16InterpolationExpHistogram
};
class TFCS2DFunctionInt8Int16Int32InterpolationExpHistogram
    : public TFCS2DFunctionTemplateInterpolationExpHistogram<
          TFCS1DFunction_HistogramInt8BinEdges,
          TFCS1DFunction_HistogramInt16BinEdges,
          uint32_t,
          float>
{
public:
  TFCS2DFunctionInt8Int16Int32InterpolationExpHistogram(TH2* h = nullptr)
      : TFCS2DFunctionTemplateInterpolationExpHistogram<
            TFCS1DFunction_HistogramInt8BinEdges,
            TFCS1DFunction_HistogramInt16BinEdges,
            uint32_t,
            float>(h) {};

  ClassDef(TFCS2DFunctionInt8Int16Int32InterpolationExpHistogram,
           1)  // TFCS2DFunctionInt8Int16Int32InterpolationExpHistogram
};

class TFCS2DFunctionInt8Int32Int8InterpolationExpHistogram
    : public TFCS2DFunctionTemplateInterpolationExpHistogram<
          TFCS1DFunction_HistogramInt8BinEdges,
          TFCS1DFunction_HistogramInt32BinEdges,
          uint8_t,
          float>
{
public:
  TFCS2DFunctionInt8Int32Int8InterpolationExpHistogram(TH2* h = nullptr)
      : TFCS2DFunctionTemplateInterpolationExpHistogram<
            TFCS1DFunction_HistogramInt8BinEdges,
            TFCS1DFunction_HistogramInt32BinEdges,
            uint8_t,
            float>(h) {};

  ClassDef(TFCS2DFunctionInt8Int32Int8InterpolationExpHistogram,
           1)  // TFCS2DFunctionInt8Int32Int8InterpolationExpHistogram
};
class TFCS2DFunctionInt8Int32Int16InterpolationExpHistogram
    : public TFCS2DFunctionTemplateInterpolationExpHistogram<
          TFCS1DFunction_HistogramInt8BinEdges,
          TFCS1DFunction_HistogramInt32BinEdges,
          uint16_t,
          float>
{
public:
  TFCS2DFunctionInt8Int32Int16InterpolationExpHistogram(TH2* h = nullptr)
      : TFCS2DFunctionTemplateInterpolationExpHistogram<
            TFCS1DFunction_HistogramInt8BinEdges,
            TFCS1DFunction_HistogramInt32BinEdges,
            uint16_t,
            float>(h) {};

  ClassDef(TFCS2DFunctionInt8Int32Int16InterpolationExpHistogram,
           1)  // TFCS2DFunctionInt8Int32Int16InterpolationExpHistogram
};
class TFCS2DFunctionInt8Int32Int32InterpolationExpHistogram
    : public TFCS2DFunctionTemplateInterpolationExpHistogram<
          TFCS1DFunction_HistogramInt8BinEdges,
          TFCS1DFunction_HistogramInt32BinEdges,
          uint32_t,
          float>
{
public:
  TFCS2DFunctionInt8Int32Int32InterpolationExpHistogram(TH2* h = nullptr)
      : TFCS2DFunctionTemplateInterpolationExpHistogram<
            TFCS1DFunction_HistogramInt8BinEdges,
            TFCS1DFunction_HistogramInt32BinEdges,
            uint32_t,
            float>(h) {};

  ClassDef(TFCS2DFunctionInt8Int32Int32InterpolationExpHistogram,
           1)  // TFCS2DFunctionInt8Int32Int32InterpolationExpHistogram
};

class TFCS2DFunctionInt16Int8Int8InterpolationExpHistogram
    : public TFCS2DFunctionTemplateInterpolationExpHistogram<
          TFCS1DFunction_HistogramInt16BinEdges,
          TFCS1DFunction_HistogramInt8BinEdges,
          uint8_t,
          float>
{
public:
  TFCS2DFunctionInt16Int8Int8InterpolationExpHistogram(TH2* h = nullptr)
      : TFCS2DFunctionTemplateInterpolationExpHistogram<
            TFCS1DFunction_HistogramInt16BinEdges,
            TFCS1DFunction_HistogramInt8BinEdges,
            uint8_t,
            float>(h) {};

  ClassDef(TFCS2DFunctionInt16Int8Int8InterpolationExpHistogram,
           1)  // TFCS2DFunctionInt16Int8Int8InterpolationExpHistogram
};
class TFCS2DFunctionInt16Int8Int16InterpolationExpHistogram
    : public TFCS2DFunctionTemplateInterpolationExpHistogram<
          TFCS1DFunction_HistogramInt16BinEdges,
          TFCS1DFunction_HistogramInt8BinEdges,
          uint16_t,
          float>
{
public:
  TFCS2DFunctionInt16Int8Int16InterpolationExpHistogram(TH2* h = nullptr)
      : TFCS2DFunctionTemplateInterpolationExpHistogram<
            TFCS1DFunction_HistogramInt16BinEdges,
            TFCS1DFunction_HistogramInt8BinEdges,
            uint16_t,
            float>(h) {};

  ClassDef(TFCS2DFunctionInt16Int8Int16InterpolationExpHistogram,
           1)  // TFCS2DFunctionInt16Int8Int16InterpolationExpHistogram
};
class TFCS2DFunctionInt16Int8Int32InterpolationExpHistogram
    : public TFCS2DFunctionTemplateInterpolationExpHistogram<
          TFCS1DFunction_HistogramInt16BinEdges,
          TFCS1DFunction_HistogramInt8BinEdges,
          uint32_t,
          float>
{
public:
  TFCS2DFunctionInt16Int8Int32InterpolationExpHistogram(TH2* h = nullptr)
      : TFCS2DFunctionTemplateInterpolationExpHistogram<
            TFCS1DFunction_HistogramInt16BinEdges,
            TFCS1DFunction_HistogramInt8BinEdges,
            uint32_t,
            float>(h) {};

  ClassDef(TFCS2DFunctionInt16Int8Int32InterpolationExpHistogram,
           1)  // TFCS2DFunctionInt16Int8Int32InterpolationExpHistogram
};

////

class TFCS2DFunctionInt16Int16Int8InterpolationExpHistogram
    : public TFCS2DFunctionTemplateInterpolationExpHistogram<
          TFCS1DFunction_HistogramInt16BinEdges,
          TFCS1DFunction_HistogramInt16BinEdges,
          uint8_t,
          float>
{
public:
  TFCS2DFunctionInt16Int16Int8InterpolationExpHistogram(TH2* h = nullptr)
      : TFCS2DFunctionTemplateInterpolationExpHistogram<
            TFCS1DFunction_HistogramInt16BinEdges,
            TFCS1DFunction_HistogramInt16BinEdges,
            uint8_t,
            float>(h) {};

  ClassDef(TFCS2DFunctionInt16Int16Int8InterpolationExpHistogram,
           1)  // TFCS2DFunctionInt16Int16Int8InterpolationExpHistogram
};
class TFCS2DFunctionInt16Int16Int16InterpolationExpHistogram
    : public TFCS2DFunctionTemplateInterpolationExpHistogram<
          TFCS1DFunction_HistogramInt16BinEdges,
          TFCS1DFunction_HistogramInt16BinEdges,
          uint16_t,
          float>
{
public:
  TFCS2DFunctionInt16Int16Int16InterpolationExpHistogram(TH2* h = nullptr)
      : TFCS2DFunctionTemplateInterpolationExpHistogram<
            TFCS1DFunction_HistogramInt16BinEdges,
            TFCS1DFunction_HistogramInt16BinEdges,
            uint16_t,
            float>(h) {};

  ClassDef(TFCS2DFunctionInt16Int16Int16InterpolationExpHistogram,
           1)  // TFCS2DFunctionInt16Int16Int16InterpolationExpHistogram
};
class TFCS2DFunctionInt16Int16Int32InterpolationExpHistogram
    : public TFCS2DFunctionTemplateInterpolationExpHistogram<
          TFCS1DFunction_HistogramInt16BinEdges,
          TFCS1DFunction_HistogramInt16BinEdges,
          uint32_t,
          float>
{
public:
  TFCS2DFunctionInt16Int16Int32InterpolationExpHistogram(TH2* h = nullptr)
      : TFCS2DFunctionTemplateInterpolationExpHistogram<
            TFCS1DFunction_HistogramInt16BinEdges,
            TFCS1DFunction_HistogramInt16BinEdges,
            uint32_t,
            float>(h) {};

  ClassDef(TFCS2DFunctionInt16Int16Int32InterpolationExpHistogram,
           1)  // TFCS2DFunctionInt16Int16Int32InterpolationExpHistogram
};

class TFCS2DFunctionInt16Int32Int8InterpolationExpHistogram
    : public TFCS2DFunctionTemplateInterpolationExpHistogram<
          TFCS1DFunction_HistogramInt16BinEdges,
          TFCS1DFunction_HistogramInt32BinEdges,
          uint8_t,
          float>
{
public:
  TFCS2DFunctionInt16Int32Int8InterpolationExpHistogram(TH2* h = nullptr)
      : TFCS2DFunctionTemplateInterpolationExpHistogram<
            TFCS1DFunction_HistogramInt16BinEdges,
            TFCS1DFunction_HistogramInt32BinEdges,
            uint8_t,
            float>(h) {};

  ClassDef(TFCS2DFunctionInt16Int32Int8InterpolationExpHistogram,
           1)  // TFCS2DFunctionInt16Int32Int8InterpolationExpHistogram
};
class TFCS2DFunctionInt16Int32Int16InterpolationExpHistogram
    : public TFCS2DFunctionTemplateInterpolationExpHistogram<
          TFCS1DFunction_HistogramInt16BinEdges,
          TFCS1DFunction_HistogramInt32BinEdges,
          uint16_t,
          float>
{
public:
  TFCS2DFunctionInt16Int32Int16InterpolationExpHistogram(TH2* h = nullptr)
      : TFCS2DFunctionTemplateInterpolationExpHistogram<
            TFCS1DFunction_HistogramInt16BinEdges,
            TFCS1DFunction_HistogramInt32BinEdges,
            uint16_t,
            float>(h) {};

  ClassDef(TFCS2DFunctionInt16Int32Int16InterpolationExpHistogram,
           1)  // TFCS2DFunctionInt16Int32Int16InterpolationExpHistogram
};
class TFCS2DFunctionInt16Int32Int32InterpolationExpHistogram
    : public TFCS2DFunctionTemplateInterpolationExpHistogram<
          TFCS1DFunction_HistogramInt16BinEdges,
          TFCS1DFunction_HistogramInt32BinEdges,
          uint32_t,
          float>
{
public:
  TFCS2DFunctionInt16Int32Int32InterpolationExpHistogram(TH2* h = nullptr)
      : TFCS2DFunctionTemplateInterpolationExpHistogram<
            TFCS1DFunction_HistogramInt16BinEdges,
            TFCS1DFunction_HistogramInt32BinEdges,
            uint32_t,
            float>(h) {};

  ClassDef(TFCS2DFunctionInt16Int32Int32InterpolationExpHistogram,
           1)  // TFCS2DFunctionInt16Int32Int32InterpolationExpHistogram
};

class TFCS2DFunctionInt32Int8Int8InterpolationExpHistogram
    : public TFCS2DFunctionTemplateInterpolationExpHistogram<
          TFCS1DFunction_HistogramInt32BinEdges,
          TFCS1DFunction_HistogramInt8BinEdges,
          uint8_t,
          float>
{
public:
  TFCS2DFunctionInt32Int8Int8InterpolationExpHistogram(TH2* h = nullptr)
      : TFCS2DFunctionTemplateInterpolationExpHistogram<
            TFCS1DFunction_HistogramInt32BinEdges,
            TFCS1DFunction_HistogramInt8BinEdges,
            uint8_t,
            float>(h) {};

  ClassDef(TFCS2DFunctionInt32Int8Int8InterpolationExpHistogram,
           1)  // TFCS2DFunctionInt32Int8Int8InterpolationExpHistogram
};
class TFCS2DFunctionInt32Int8Int16InterpolationExpHistogram
    : public TFCS2DFunctionTemplateInterpolationExpHistogram<
          TFCS1DFunction_HistogramInt32BinEdges,
          TFCS1DFunction_HistogramInt8BinEdges,
          uint16_t,
          float>
{
public:
  TFCS2DFunctionInt32Int8Int16InterpolationExpHistogram(TH2* h = nullptr)
      : TFCS2DFunctionTemplateInterpolationExpHistogram<
            TFCS1DFunction_HistogramInt32BinEdges,
            TFCS1DFunction_HistogramInt8BinEdges,
            uint16_t,
            float>(h) {};

  ClassDef(TFCS2DFunctionInt32Int8Int16InterpolationExpHistogram,
           1)  // TFCS2DFunctionInt32Int8Int16InterpolationExpHistogram
};
class TFCS2DFunctionInt32Int8Int32InterpolationExpHistogram
    : public TFCS2DFunctionTemplateInterpolationExpHistogram<
          TFCS1DFunction_HistogramInt32BinEdges,
          TFCS1DFunction_HistogramInt8BinEdges,
          uint32_t,
          float>
{
public:
  TFCS2DFunctionInt32Int8Int32InterpolationExpHistogram(TH2* h = nullptr)
      : TFCS2DFunctionTemplateInterpolationExpHistogram<
            TFCS1DFunction_HistogramInt32BinEdges,
            TFCS1DFunction_HistogramInt8BinEdges,
            uint32_t,
            float>(h) {};

  ClassDef(TFCS2DFunctionInt32Int8Int32InterpolationExpHistogram,
           1)  // TFCS2DFunctionInt32Int8Int32InterpolationExpHistogram
};

////

class TFCS2DFunctionInt32Int16Int8InterpolationExpHistogram
    : public TFCS2DFunctionTemplateInterpolationExpHistogram<
          TFCS1DFunction_HistogramInt32BinEdges,
          TFCS1DFunction_HistogramInt16BinEdges,
          uint8_t,
          float>
{
public:
  TFCS2DFunctionInt32Int16Int8InterpolationExpHistogram(TH2* h = nullptr)
      : TFCS2DFunctionTemplateInterpolationExpHistogram<
            TFCS1DFunction_HistogramInt32BinEdges,
            TFCS1DFunction_HistogramInt16BinEdges,
            uint8_t,
            float>(h) {};

  ClassDef(TFCS2DFunctionInt32Int16Int8InterpolationExpHistogram,
           1)  // TFCS2DFunctionInt32Int16Int8InterpolationExpHistogram
};
class TFCS2DFunctionInt32Int16Int16InterpolationExpHistogram
    : public TFCS2DFunctionTemplateInterpolationExpHistogram<
          TFCS1DFunction_HistogramInt32BinEdges,
          TFCS1DFunction_HistogramInt16BinEdges,
          uint16_t,
          float>
{
public:
  TFCS2DFunctionInt32Int16Int16InterpolationExpHistogram(TH2* h = nullptr)
      : TFCS2DFunctionTemplateInterpolationExpHistogram<
            TFCS1DFunction_HistogramInt32BinEdges,
            TFCS1DFunction_HistogramInt16BinEdges,
            uint16_t,
            float>(h) {};

  ClassDef(TFCS2DFunctionInt32Int16Int16InterpolationExpHistogram,
           1)  // TFCS2DFunctionInt32Int16Int16InterpolationExpHistogram
};
class TFCS2DFunctionInt32Int16Int32InterpolationExpHistogram
    : public TFCS2DFunctionTemplateInterpolationExpHistogram<
          TFCS1DFunction_HistogramInt32BinEdges,
          TFCS1DFunction_HistogramInt16BinEdges,
          uint32_t,
          float>
{
public:
  TFCS2DFunctionInt32Int16Int32InterpolationExpHistogram(TH2* h = nullptr)
      : TFCS2DFunctionTemplateInterpolationExpHistogram<
            TFCS1DFunction_HistogramInt32BinEdges,
            TFCS1DFunction_HistogramInt16BinEdges,
            uint32_t,
            float>(h) {};

  ClassDef(TFCS2DFunctionInt32Int16Int32InterpolationExpHistogram,
           1)  // TFCS2DFunctionInt32Int16Int32InterpolationExpHistogram
};

class TFCS2DFunctionInt32Int32Int8InterpolationExpHistogram
    : public TFCS2DFunctionTemplateInterpolationExpHistogram<
          TFCS1DFunction_HistogramInt32BinEdges,
          TFCS1DFunction_HistogramInt32BinEdges,
          uint8_t,
          float>
{
public:
  TFCS2DFunctionInt32Int32Int8InterpolationExpHistogram(TH2* h = nullptr)
      : TFCS2DFunctionTemplateInterpolationExpHistogram<
            TFCS1DFunction_HistogramInt32BinEdges,
            TFCS1DFunction_HistogramInt32BinEdges,
            uint8_t,
            float>(h) {};

  ClassDef(TFCS2DFunctionInt32Int32Int8InterpolationExpHistogram,
           1)  // TFCS2DFunctionInt32Int32Int8InterpolationExpHistogram
};
class TFCS2DFunctionInt32Int32Int16InterpolationExpHistogram
    : public TFCS2DFunctionTemplateInterpolationExpHistogram<
          TFCS1DFunction_HistogramInt32BinEdges,
          TFCS1DFunction_HistogramInt32BinEdges,
          uint16_t,
          float>
{
public:
  TFCS2DFunctionInt32Int32Int16InterpolationExpHistogram(TH2* h = nullptr)
      : TFCS2DFunctionTemplateInterpolationExpHistogram<
            TFCS1DFunction_HistogramInt32BinEdges,
            TFCS1DFunction_HistogramInt32BinEdges,
            uint16_t,
            float>(h) {};

  ClassDef(TFCS2DFunctionInt32Int32Int16InterpolationExpHistogram,
           1)  // TFCS2DFunctionInt32Int32Int16InterpolationExpHistogram
};
class TFCS2DFunctionInt32Int32Int32InterpolationExpHistogram
    : public TFCS2DFunctionTemplateInterpolationExpHistogram<
          TFCS1DFunction_HistogramInt32BinEdges,
          TFCS1DFunction_HistogramInt32BinEdges,
          uint32_t,
          float>
{
public:
  TFCS2DFunctionInt32Int32Int32InterpolationExpHistogram(TH2* h = nullptr)
      : TFCS2DFunctionTemplateInterpolationExpHistogram<
            TFCS1DFunction_HistogramInt32BinEdges,
            TFCS1DFunction_HistogramInt32BinEdges,
            uint32_t,
            float>(h) {};

  ClassDef(TFCS2DFunctionInt32Int32Int32InterpolationExpHistogram,
           1)  // TFCS2DFunctionInt32Int32Int32InterpolationExpHistogram
};

}  // namespace FastCaloSim::Core

// clang-format off
#if defined(__ROOTCLING__) && defined(__FastCaloSimStandAlone__)
#  pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateInterpolationExpHistogram \
          < FastCaloSim::Core::TFCS1DFunction_HistogramInt8BinEdges, \
      FastCaloSim::Core::TFCS1DFunction_HistogramInt8BinEdges, uint8_t, float> + ;
#  pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateInterpolationExpHistogram \
          < FastCaloSim::Core::TFCS1DFunction_HistogramInt8BinEdges, \
      FastCaloSim::Core::TFCS1DFunction_HistogramInt8BinEdges, uint16_t, float> + ;
#  pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateInterpolationExpHistogram \
          < FastCaloSim::Core::TFCS1DFunction_HistogramInt8BinEdges, \
      FastCaloSim::Core::TFCS1DFunction_HistogramInt8BinEdges, uint32_t, float> + ;
#  pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateInterpolationExpHistogram \
          < FastCaloSim::Core::TFCS1DFunction_HistogramInt8BinEdges, \
      FastCaloSim::Core::TFCS1DFunction_HistogramInt16BinEdges, uint8_t, float> + ;
#  pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateInterpolationExpHistogram \
          < FastCaloSim::Core::TFCS1DFunction_HistogramInt8BinEdges, \
      FastCaloSim::Core::TFCS1DFunction_HistogramInt16BinEdges, uint16_t, float> + ;
#  pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateInterpolationExpHistogram \
          < FastCaloSim::Core::TFCS1DFunction_HistogramInt8BinEdges, \
      FastCaloSim::Core::TFCS1DFunction_HistogramInt16BinEdges, uint32_t, float> + ;
#  pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateInterpolationExpHistogram \
          < FastCaloSim::Core::TFCS1DFunction_HistogramInt8BinEdges, \
      FastCaloSim::Core::TFCS1DFunction_HistogramInt32BinEdges, uint8_t, float> + ;
#  pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateInterpolationExpHistogram \
          < FastCaloSim::Core::TFCS1DFunction_HistogramInt8BinEdges, \
      FastCaloSim::Core::TFCS1DFunction_HistogramInt32BinEdges, uint16_t, float> + ;
#  pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateInterpolationExpHistogram \
          < FastCaloSim::Core::TFCS1DFunction_HistogramInt8BinEdges, \
      FastCaloSim::Core::TFCS1DFunction_HistogramInt32BinEdges, uint32_t, float> + ;
#  pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateInterpolationExpHistogram \
          < FastCaloSim::Core::TFCS1DFunction_HistogramInt16BinEdges, \
      FastCaloSim::Core::TFCS1DFunction_HistogramInt8BinEdges, uint8_t, float> + ;
#  pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateInterpolationExpHistogram \
          < FastCaloSim::Core::TFCS1DFunction_HistogramInt16BinEdges, \
      FastCaloSim::Core::TFCS1DFunction_HistogramInt8BinEdges, uint16_t, float> + ;
#  pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateInterpolationExpHistogram \
          < FastCaloSim::Core::TFCS1DFunction_HistogramInt16BinEdges, \
      FastCaloSim::Core::TFCS1DFunction_HistogramInt8BinEdges, uint32_t, float> + ;
#  pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateInterpolationExpHistogram \
          < FastCaloSim::Core::TFCS1DFunction_HistogramInt16BinEdges, \
      FastCaloSim::Core::TFCS1DFunction_HistogramInt16BinEdges, uint8_t, float> + ;
#  pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateInterpolationExpHistogram \
          < FastCaloSim::Core::TFCS1DFunction_HistogramInt16BinEdges, \
      FastCaloSim::Core::TFCS1DFunction_HistogramInt16BinEdges, uint16_t, float> + ;
#  pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateInterpolationExpHistogram \
          < FastCaloSim::Core::TFCS1DFunction_HistogramInt16BinEdges, \
      FastCaloSim::Core::TFCS1DFunction_HistogramInt16BinEdges, uint32_t, float> + ;
#  pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateInterpolationExpHistogram \
          < FastCaloSim::Core::TFCS1DFunction_HistogramInt16BinEdges, \
      FastCaloSim::Core::TFCS1DFunction_HistogramInt32BinEdges, uint8_t, float> + ;
#  pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateInterpolationExpHistogram \
          < FastCaloSim::Core::TFCS1DFunction_HistogramInt16BinEdges, \
      FastCaloSim::Core::TFCS1DFunction_HistogramInt32BinEdges, uint16_t, float> + ;
#  pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateInterpolationExpHistogram \
          < FastCaloSim::Core::TFCS1DFunction_HistogramInt16BinEdges, \
      FastCaloSim::Core::TFCS1DFunction_HistogramInt32BinEdges, uint32_t, float> + ;
#  pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateInterpolationExpHistogram \
          < FastCaloSim::Core::TFCS1DFunction_HistogramInt32BinEdges, \
      FastCaloSim::Core::TFCS1DFunction_HistogramInt8BinEdges, uint8_t, float> + ;
#  pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateInterpolationExpHistogram \
          < FastCaloSim::Core::TFCS1DFunction_HistogramInt32BinEdges, \
      FastCaloSim::Core::TFCS1DFunction_HistogramInt8BinEdges, uint16_t, float> + ;
#  pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateInterpolationExpHistogram \
          < FastCaloSim::Core::TFCS1DFunction_HistogramInt32BinEdges, \
      FastCaloSim::Core::TFCS1DFunction_HistogramInt8BinEdges, uint32_t, float> + ;
#  pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateInterpolationExpHistogram \
          < FastCaloSim::Core::TFCS1DFunction_HistogramInt32BinEdges, \
      FastCaloSim::Core::TFCS1DFunction_HistogramInt16BinEdges, uint8_t, float> + ;
#  pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateInterpolationExpHistogram \
          < FastCaloSim::Core::TFCS1DFunction_HistogramInt32BinEdges, \
      FastCaloSim::Core::TFCS1DFunction_HistogramInt16BinEdges, uint16_t, float> + ;
#  pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateInterpolationExpHistogram \
          < FastCaloSim::Core::TFCS1DFunction_HistogramInt32BinEdges, \
      FastCaloSim::Core::TFCS1DFunction_HistogramInt16BinEdges, uint32_t, float> + ;
#  pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateInterpolationExpHistogram \
          < FastCaloSim::Core::TFCS1DFunction_HistogramInt32BinEdges, \
      FastCaloSim::Core::TFCS1DFunction_HistogramInt32BinEdges, uint8_t, float> + ;
#  pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateInterpolationExpHistogram \
          < FastCaloSim::Core::TFCS1DFunction_HistogramInt32BinEdges, \
      FastCaloSim::Core::TFCS1DFunction_HistogramInt32BinEdges, uint16_t, float> + ;
#  pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateInterpolationExpHistogram \
          < FastCaloSim::Core::TFCS1DFunction_HistogramInt32BinEdges, \
      FastCaloSim::Core::TFCS1DFunction_HistogramInt32BinEdges, uint32_t, float> + ;

#  pragma link C++ class FastCaloSim::Core::TFCS2DFunctionInt8Int8Int8InterpolationExpHistogram + ;
#  pragma link C++ class FastCaloSim::Core::TFCS2DFunctionInt8Int8Int16InterpolationExpHistogram \
          + ;
#  pragma link C++ class FastCaloSim::Core::TFCS2DFunctionInt8Int8Int32InterpolationExpHistogram \
          + ;
#  pragma link C++ class FastCaloSim::Core::TFCS2DFunctionInt8Int16Int8InterpolationExpHistogram \
          + ;
#  pragma link C++ class FastCaloSim::Core::TFCS2DFunctionInt8Int16Int16InterpolationExpHistogram \
          + ;
#  pragma link C++ class FastCaloSim::Core::TFCS2DFunctionInt8Int16Int32InterpolationExpHistogram \
          + ;
#  pragma link C++ class FastCaloSim::Core::TFCS2DFunctionInt8Int32Int8InterpolationExpHistogram \
          + ;
#  pragma link C++ class FastCaloSim::Core::TFCS2DFunctionInt8Int32Int16InterpolationExpHistogram \
          + ;
#  pragma link C++ class FastCaloSim::Core::TFCS2DFunctionInt8Int32Int32InterpolationExpHistogram \
          + ;
#  pragma link C++ class FastCaloSim::Core::TFCS2DFunctionInt16Int8Int8InterpolationExpHistogram \
          + ;
#  pragma link C++ class FastCaloSim::Core::TFCS2DFunctionInt16Int8Int16InterpolationExpHistogram \
          + ;
#  pragma link C++ class FastCaloSim::Core::TFCS2DFunctionInt16Int8Int32InterpolationExpHistogram \
          + ;
#  pragma link C++ class FastCaloSim::Core::TFCS2DFunctionInt16Int16Int8InterpolationExpHistogram \
          + ;
#  pragma link \
          C++ class FastCaloSim::Core::TFCS2DFunctionInt16Int16Int16InterpolationExpHistogram \
          + ;
#  pragma link \
          C++ class FastCaloSim::Core::TFCS2DFunctionInt16Int16Int32InterpolationExpHistogram \
          + ;
#  pragma link C++ class FastCaloSim::Core::TFCS2DFunctionInt16Int32Int8InterpolationExpHistogram \
          + ;
#  pragma link \
          C++ class FastCaloSim::Core::TFCS2DFunctionInt16Int32Int16InterpolationExpHistogram \
          + ;
#  pragma link \
          C++ class FastCaloSim::Core::TFCS2DFunctionInt16Int32Int32InterpolationExpHistogram \
          + ;
#  pragma link C++ class FastCaloSim::Core::TFCS2DFunctionInt32Int8Int8InterpolationExpHistogram \
          + ;
#  pragma link C++ class FastCaloSim::Core::TFCS2DFunctionInt32Int8Int16InterpolationExpHistogram \
          + ;
#  pragma link C++ class FastCaloSim::Core::TFCS2DFunctionInt32Int8Int32InterpolationExpHistogram \
          + ;
#  pragma link C++ class FastCaloSim::Core::TFCS2DFunctionInt32Int16Int8InterpolationExpHistogram \
          + ;
#  pragma link \
          C++ class FastCaloSim::Core::TFCS2DFunctionInt32Int16Int16InterpolationExpHistogram \
          + ;
#  pragma link \
          C++ class FastCaloSim::Core::TFCS2DFunctionInt32Int16Int32InterpolationExpHistogram \
          + ;
#  pragma link C++ class FastCaloSim::Core::TFCS2DFunctionInt32Int32Int8InterpolationExpHistogram \
          + ;
#  pragma link \
          C++ class FastCaloSim::Core::TFCS2DFunctionInt32Int32Int16InterpolationExpHistogram \
          + ;
#  pragma link \
          C++ class FastCaloSim::Core::TFCS2DFunctionInt32Int32Int32InterpolationExpHistogram \
          + ;

#endif
// clang-format on
#endif
