// Copyright (c) 2024 CERN for the benefit of the FastCaloSim project

#ifndef TFCS2DFunctionTemplateInterpolationHistogram_h
#define TFCS2DFunctionTemplateInterpolationHistogram_h

#include <iostream>

#include "FastCaloSim/Core/TFCS2DFunctionTemplateHistogram.h"
#include "TBuffer.h"
#include "TFCS2DFunction.h"
#include "TH1.h"
#include "TH2.h"

namespace FastCaloSim::Core
{
template<typename Txvec, typename Tyvec, typename Tz, typename Trandom = float>
class TFCS2DFunctionTemplateInterpolationHistogram
    : public TFCS2DFunctionTemplateHistogram<Txvec, Tyvec, Tz, Trandom>
{
public:
  TFCS2DFunctionTemplateInterpolationHistogram(TH2* hist = nullptr)
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

  /// Function gets random number rnd in the range [0,1) as argument
  /// and returns function value according to a histogram distribution.
  /// A linear interpolation is done within each bin
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

    Trandom dfrac;
    auto l = m_HistoBordersy.GetBinLength(ibiny);
    Trandom dfracprev;
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

    Trandom dfracnext;
    auto lnext = l;
    if (ibiny < nbinsy - 1) {
      dfracnext = m_HistoContents.get_fraction(ibin + 1) - frac;
      lnext = m_HistoBordersy.GetBinLength(ibiny + 1);
    } else {
      dfracnext = dfrac;
    }

    Trandom m;
    if (dfrac > 0) {
      Trandom mnext = dfracnext / lnext - dfrac / l;
      mnext /= l + lnext;
      Trandom mprev = dfrac / l - dfracprev / lprev;
      mprev /= l + lprev;
      m = (mnext + mprev) * l * l / dfrac;
    } else
      m = 0;

    valuey = m_HistoBordersy.position_lin(ibiny, m, residual_rnd);
    valuex = m_HistoBordersx.position(ibinx, rnd2);
  }

  ClassDef(TFCS2DFunctionTemplateInterpolationHistogram,
           1)  // TFCS1DFunctionTemplateInterpolationHistogram
};

class TFCS2DFunctionInt8Int8Int8InterpolationHistogram
    : public TFCS2DFunctionTemplateInterpolationHistogram<
          TFCS1DFunction_HistogramInt8BinEdges,
          TFCS1DFunction_HistogramInt8BinEdges,
          uint8_t,
          float>
{
public:
  TFCS2DFunctionInt8Int8Int8InterpolationHistogram(TH2* h = nullptr)
      : TFCS2DFunctionTemplateInterpolationHistogram<
            TFCS1DFunction_HistogramInt8BinEdges,
            TFCS1DFunction_HistogramInt8BinEdges,
            uint8_t,
            float>(h) {};

  ClassDef(TFCS2DFunctionInt8Int8Int8InterpolationHistogram,
           1)  // TFCS2DFunctionInt8Int8Int8InterpolationHistogram
};
class TFCS2DFunctionInt8Int8Int16InterpolationHistogram
    : public TFCS2DFunctionTemplateInterpolationHistogram<
          TFCS1DFunction_HistogramInt8BinEdges,
          TFCS1DFunction_HistogramInt8BinEdges,
          uint16_t,
          float>
{
public:
  TFCS2DFunctionInt8Int8Int16InterpolationHistogram(TH2* h = nullptr)
      : TFCS2DFunctionTemplateInterpolationHistogram<
            TFCS1DFunction_HistogramInt8BinEdges,
            TFCS1DFunction_HistogramInt8BinEdges,
            uint16_t,
            float>(h) {};

  ClassDef(TFCS2DFunctionInt8Int8Int16InterpolationHistogram,
           1)  // TFCS2DFunctionInt8Int8Int16InterpolationHistogram
};
class TFCS2DFunctionInt8Int8Int32InterpolationHistogram
    : public TFCS2DFunctionTemplateInterpolationHistogram<
          TFCS1DFunction_HistogramInt8BinEdges,
          TFCS1DFunction_HistogramInt8BinEdges,
          uint32_t,
          float>
{
public:
  TFCS2DFunctionInt8Int8Int32InterpolationHistogram(TH2* h = nullptr)
      : TFCS2DFunctionTemplateInterpolationHistogram<
            TFCS1DFunction_HistogramInt8BinEdges,
            TFCS1DFunction_HistogramInt8BinEdges,
            uint32_t,
            float>(h) {};

  ClassDef(TFCS2DFunctionInt8Int8Int32InterpolationHistogram,
           1)  // TFCS2DFunctionInt8Int8Int32InterpolationHistogram
};

////

class TFCS2DFunctionInt8Int16Int8InterpolationHistogram
    : public TFCS2DFunctionTemplateInterpolationHistogram<
          TFCS1DFunction_HistogramInt8BinEdges,
          TFCS1DFunction_HistogramInt16BinEdges,
          uint8_t,
          float>
{
public:
  TFCS2DFunctionInt8Int16Int8InterpolationHistogram(TH2* h = nullptr)
      : TFCS2DFunctionTemplateInterpolationHistogram<
            TFCS1DFunction_HistogramInt8BinEdges,
            TFCS1DFunction_HistogramInt16BinEdges,
            uint8_t,
            float>(h) {};

  ClassDef(TFCS2DFunctionInt8Int16Int8InterpolationHistogram,
           1)  // TFCS2DFunctionInt8Int16Int8InterpolationHistogram
};
class TFCS2DFunctionInt8Int16Int16InterpolationHistogram
    : public TFCS2DFunctionTemplateInterpolationHistogram<
          TFCS1DFunction_HistogramInt8BinEdges,
          TFCS1DFunction_HistogramInt16BinEdges,
          uint16_t,
          float>
{
public:
  TFCS2DFunctionInt8Int16Int16InterpolationHistogram(TH2* h = nullptr)
      : TFCS2DFunctionTemplateInterpolationHistogram<
            TFCS1DFunction_HistogramInt8BinEdges,
            TFCS1DFunction_HistogramInt16BinEdges,
            uint16_t,
            float>(h) {};

  ClassDef(TFCS2DFunctionInt8Int16Int16InterpolationHistogram,
           1)  // TFCS2DFunctionInt8Int16Int16InterpolationHistogram
};
class TFCS2DFunctionInt8Int16Int32InterpolationHistogram
    : public TFCS2DFunctionTemplateInterpolationHistogram<
          TFCS1DFunction_HistogramInt8BinEdges,
          TFCS1DFunction_HistogramInt16BinEdges,
          uint32_t,
          float>
{
public:
  TFCS2DFunctionInt8Int16Int32InterpolationHistogram(TH2* h = nullptr)
      : TFCS2DFunctionTemplateInterpolationHistogram<
            TFCS1DFunction_HistogramInt8BinEdges,
            TFCS1DFunction_HistogramInt16BinEdges,
            uint32_t,
            float>(h) {};

  ClassDef(TFCS2DFunctionInt8Int16Int32InterpolationHistogram,
           1)  // TFCS2DFunctionInt8Int16Int32InterpolationHistogram
};

class TFCS2DFunctionInt8Int32Int8InterpolationHistogram
    : public TFCS2DFunctionTemplateInterpolationHistogram<
          TFCS1DFunction_HistogramInt8BinEdges,
          TFCS1DFunction_HistogramInt32BinEdges,
          uint8_t,
          float>
{
public:
  TFCS2DFunctionInt8Int32Int8InterpolationHistogram(TH2* h = nullptr)
      : TFCS2DFunctionTemplateInterpolationHistogram<
            TFCS1DFunction_HistogramInt8BinEdges,
            TFCS1DFunction_HistogramInt32BinEdges,
            uint8_t,
            float>(h) {};

  ClassDef(TFCS2DFunctionInt8Int32Int8InterpolationHistogram,
           1)  // TFCS2DFunctionInt8Int32Int8InterpolationHistogram
};
class TFCS2DFunctionInt8Int32Int16InterpolationHistogram
    : public TFCS2DFunctionTemplateInterpolationHistogram<
          TFCS1DFunction_HistogramInt8BinEdges,
          TFCS1DFunction_HistogramInt32BinEdges,
          uint16_t,
          float>
{
public:
  TFCS2DFunctionInt8Int32Int16InterpolationHistogram(TH2* h = nullptr)
      : TFCS2DFunctionTemplateInterpolationHistogram<
            TFCS1DFunction_HistogramInt8BinEdges,
            TFCS1DFunction_HistogramInt32BinEdges,
            uint16_t,
            float>(h) {};

  ClassDef(TFCS2DFunctionInt8Int32Int16InterpolationHistogram,
           1)  // TFCS2DFunctionInt8Int32Int16InterpolationHistogram
};
class TFCS2DFunctionInt8Int32Int32InterpolationHistogram
    : public TFCS2DFunctionTemplateInterpolationHistogram<
          TFCS1DFunction_HistogramInt8BinEdges,
          TFCS1DFunction_HistogramInt32BinEdges,
          uint32_t,
          float>
{
public:
  TFCS2DFunctionInt8Int32Int32InterpolationHistogram(TH2* h = nullptr)
      : TFCS2DFunctionTemplateInterpolationHistogram<
            TFCS1DFunction_HistogramInt8BinEdges,
            TFCS1DFunction_HistogramInt32BinEdges,
            uint32_t,
            float>(h) {};

  ClassDef(TFCS2DFunctionInt8Int32Int32InterpolationHistogram,
           1)  // TFCS2DFunctionInt8Int32Int32InterpolationHistogram
};

class TFCS2DFunctionInt16Int8Int8InterpolationHistogram
    : public TFCS2DFunctionTemplateInterpolationHistogram<
          TFCS1DFunction_HistogramInt16BinEdges,
          TFCS1DFunction_HistogramInt8BinEdges,
          uint8_t,
          float>
{
public:
  TFCS2DFunctionInt16Int8Int8InterpolationHistogram(TH2* h = nullptr)
      : TFCS2DFunctionTemplateInterpolationHistogram<
            TFCS1DFunction_HistogramInt16BinEdges,
            TFCS1DFunction_HistogramInt8BinEdges,
            uint8_t,
            float>(h) {};

  ClassDef(TFCS2DFunctionInt16Int8Int8InterpolationHistogram,
           1)  // TFCS2DFunctionInt16Int8Int8InterpolationHistogram
};
class TFCS2DFunctionInt16Int8Int16InterpolationHistogram
    : public TFCS2DFunctionTemplateInterpolationHistogram<
          TFCS1DFunction_HistogramInt16BinEdges,
          TFCS1DFunction_HistogramInt8BinEdges,
          uint16_t,
          float>
{
public:
  TFCS2DFunctionInt16Int8Int16InterpolationHistogram(TH2* h = nullptr)
      : TFCS2DFunctionTemplateInterpolationHistogram<
            TFCS1DFunction_HistogramInt16BinEdges,
            TFCS1DFunction_HistogramInt8BinEdges,
            uint16_t,
            float>(h) {};

  ClassDef(TFCS2DFunctionInt16Int8Int16InterpolationHistogram,
           1)  // TFCS2DFunctionInt16Int8Int16InterpolationHistogram
};
class TFCS2DFunctionInt16Int8Int32InterpolationHistogram
    : public TFCS2DFunctionTemplateInterpolationHistogram<
          TFCS1DFunction_HistogramInt16BinEdges,
          TFCS1DFunction_HistogramInt8BinEdges,
          uint32_t,
          float>
{
public:
  TFCS2DFunctionInt16Int8Int32InterpolationHistogram(TH2* h = nullptr)
      : TFCS2DFunctionTemplateInterpolationHistogram<
            TFCS1DFunction_HistogramInt16BinEdges,
            TFCS1DFunction_HistogramInt8BinEdges,
            uint32_t,
            float>(h) {};

  ClassDef(TFCS2DFunctionInt16Int8Int32InterpolationHistogram,
           1)  // TFCS2DFunctionInt16Int8Int32InterpolationHistogram
};

////

class TFCS2DFunctionInt16Int16Int8InterpolationHistogram
    : public TFCS2DFunctionTemplateInterpolationHistogram<
          TFCS1DFunction_HistogramInt16BinEdges,
          TFCS1DFunction_HistogramInt16BinEdges,
          uint8_t,
          float>
{
public:
  TFCS2DFunctionInt16Int16Int8InterpolationHistogram(TH2* h = nullptr)
      : TFCS2DFunctionTemplateInterpolationHistogram<
            TFCS1DFunction_HistogramInt16BinEdges,
            TFCS1DFunction_HistogramInt16BinEdges,
            uint8_t,
            float>(h) {};

  ClassDef(TFCS2DFunctionInt16Int16Int8InterpolationHistogram,
           1)  // TFCS2DFunctionInt16Int16Int8InterpolationHistogram
};
class TFCS2DFunctionInt16Int16Int16InterpolationHistogram
    : public TFCS2DFunctionTemplateInterpolationHistogram<
          TFCS1DFunction_HistogramInt16BinEdges,
          TFCS1DFunction_HistogramInt16BinEdges,
          uint16_t,
          float>
{
public:
  TFCS2DFunctionInt16Int16Int16InterpolationHistogram(TH2* h = nullptr)
      : TFCS2DFunctionTemplateInterpolationHistogram<
            TFCS1DFunction_HistogramInt16BinEdges,
            TFCS1DFunction_HistogramInt16BinEdges,
            uint16_t,
            float>(h) {};

  ClassDef(TFCS2DFunctionInt16Int16Int16InterpolationHistogram,
           1)  // TFCS2DFunctionInt16Int16Int16InterpolationHistogram
};
class TFCS2DFunctionInt16Int16Int32InterpolationHistogram
    : public TFCS2DFunctionTemplateInterpolationHistogram<
          TFCS1DFunction_HistogramInt16BinEdges,
          TFCS1DFunction_HistogramInt16BinEdges,
          uint32_t,
          float>
{
public:
  TFCS2DFunctionInt16Int16Int32InterpolationHistogram(TH2* h = nullptr)
      : TFCS2DFunctionTemplateInterpolationHistogram<
            TFCS1DFunction_HistogramInt16BinEdges,
            TFCS1DFunction_HistogramInt16BinEdges,
            uint32_t,
            float>(h) {};

  ClassDef(TFCS2DFunctionInt16Int16Int32InterpolationHistogram,
           1)  // TFCS2DFunctionInt16Int16Int32InterpolationHistogram
};

class TFCS2DFunctionInt16Int32Int8InterpolationHistogram
    : public TFCS2DFunctionTemplateInterpolationHistogram<
          TFCS1DFunction_HistogramInt16BinEdges,
          TFCS1DFunction_HistogramInt32BinEdges,
          uint8_t,
          float>
{
public:
  TFCS2DFunctionInt16Int32Int8InterpolationHistogram(TH2* h = nullptr)
      : TFCS2DFunctionTemplateInterpolationHistogram<
            TFCS1DFunction_HistogramInt16BinEdges,
            TFCS1DFunction_HistogramInt32BinEdges,
            uint8_t,
            float>(h) {};

  ClassDef(TFCS2DFunctionInt16Int32Int8InterpolationHistogram,
           1)  // TFCS2DFunctionInt16Int32Int8InterpolationHistogram
};
class TFCS2DFunctionInt16Int32Int16InterpolationHistogram
    : public TFCS2DFunctionTemplateInterpolationHistogram<
          TFCS1DFunction_HistogramInt16BinEdges,
          TFCS1DFunction_HistogramInt32BinEdges,
          uint16_t,
          float>
{
public:
  TFCS2DFunctionInt16Int32Int16InterpolationHistogram(TH2* h = nullptr)
      : TFCS2DFunctionTemplateInterpolationHistogram<
            TFCS1DFunction_HistogramInt16BinEdges,
            TFCS1DFunction_HistogramInt32BinEdges,
            uint16_t,
            float>(h) {};

  ClassDef(TFCS2DFunctionInt16Int32Int16InterpolationHistogram,
           1)  // TFCS2DFunctionInt16Int32Int16InterpolationHistogram
};
class TFCS2DFunctionInt16Int32Int32InterpolationHistogram
    : public TFCS2DFunctionTemplateInterpolationHistogram<
          TFCS1DFunction_HistogramInt16BinEdges,
          TFCS1DFunction_HistogramInt32BinEdges,
          uint32_t,
          float>
{
public:
  TFCS2DFunctionInt16Int32Int32InterpolationHistogram(TH2* h = nullptr)
      : TFCS2DFunctionTemplateInterpolationHistogram<
            TFCS1DFunction_HistogramInt16BinEdges,
            TFCS1DFunction_HistogramInt32BinEdges,
            uint32_t,
            float>(h) {};

  ClassDef(TFCS2DFunctionInt16Int32Int32InterpolationHistogram,
           1)  // TFCS2DFunctionInt16Int32Int32InterpolationHistogram
};

class TFCS2DFunctionInt32Int8Int8InterpolationHistogram
    : public TFCS2DFunctionTemplateInterpolationHistogram<
          TFCS1DFunction_HistogramInt32BinEdges,
          TFCS1DFunction_HistogramInt8BinEdges,
          uint8_t,
          float>
{
public:
  TFCS2DFunctionInt32Int8Int8InterpolationHistogram(TH2* h = nullptr)
      : TFCS2DFunctionTemplateInterpolationHistogram<
            TFCS1DFunction_HistogramInt32BinEdges,
            TFCS1DFunction_HistogramInt8BinEdges,
            uint8_t,
            float>(h) {};

  ClassDef(TFCS2DFunctionInt32Int8Int8InterpolationHistogram,
           1)  // TFCS2DFunctionInt32Int8Int8InterpolationHistogram
};
class TFCS2DFunctionInt32Int8Int16InterpolationHistogram
    : public TFCS2DFunctionTemplateInterpolationHistogram<
          TFCS1DFunction_HistogramInt32BinEdges,
          TFCS1DFunction_HistogramInt8BinEdges,
          uint16_t,
          float>
{
public:
  TFCS2DFunctionInt32Int8Int16InterpolationHistogram(TH2* h = nullptr)
      : TFCS2DFunctionTemplateInterpolationHistogram<
            TFCS1DFunction_HistogramInt32BinEdges,
            TFCS1DFunction_HistogramInt8BinEdges,
            uint16_t,
            float>(h) {};

  ClassDef(TFCS2DFunctionInt32Int8Int16InterpolationHistogram,
           1)  // TFCS2DFunctionInt32Int8Int16InterpolationHistogram
};
class TFCS2DFunctionInt32Int8Int32InterpolationHistogram
    : public TFCS2DFunctionTemplateInterpolationHistogram<
          TFCS1DFunction_HistogramInt32BinEdges,
          TFCS1DFunction_HistogramInt8BinEdges,
          uint32_t,
          float>
{
public:
  TFCS2DFunctionInt32Int8Int32InterpolationHistogram(TH2* h = nullptr)
      : TFCS2DFunctionTemplateInterpolationHistogram<
            TFCS1DFunction_HistogramInt32BinEdges,
            TFCS1DFunction_HistogramInt8BinEdges,
            uint32_t,
            float>(h) {};

  ClassDef(TFCS2DFunctionInt32Int8Int32InterpolationHistogram,
           1)  // TFCS2DFunctionInt32Int8Int32InterpolationHistogram
};

////

class TFCS2DFunctionInt32Int16Int8InterpolationHistogram
    : public TFCS2DFunctionTemplateInterpolationHistogram<
          TFCS1DFunction_HistogramInt32BinEdges,
          TFCS1DFunction_HistogramInt16BinEdges,
          uint8_t,
          float>
{
public:
  TFCS2DFunctionInt32Int16Int8InterpolationHistogram(TH2* h = nullptr)
      : TFCS2DFunctionTemplateInterpolationHistogram<
            TFCS1DFunction_HistogramInt32BinEdges,
            TFCS1DFunction_HistogramInt16BinEdges,
            uint8_t,
            float>(h) {};

  ClassDef(TFCS2DFunctionInt32Int16Int8InterpolationHistogram,
           1)  // TFCS2DFunctionInt32Int16Int8InterpolationHistogram
};
class TFCS2DFunctionInt32Int16Int16InterpolationHistogram
    : public TFCS2DFunctionTemplateInterpolationHistogram<
          TFCS1DFunction_HistogramInt32BinEdges,
          TFCS1DFunction_HistogramInt16BinEdges,
          uint16_t,
          float>
{
public:
  TFCS2DFunctionInt32Int16Int16InterpolationHistogram(TH2* h = nullptr)
      : TFCS2DFunctionTemplateInterpolationHistogram<
            TFCS1DFunction_HistogramInt32BinEdges,
            TFCS1DFunction_HistogramInt16BinEdges,
            uint16_t,
            float>(h) {};

  ClassDef(TFCS2DFunctionInt32Int16Int16InterpolationHistogram,
           1)  // TFCS2DFunctionInt32Int16Int16InterpolationHistogram
};
class TFCS2DFunctionInt32Int16Int32InterpolationHistogram
    : public TFCS2DFunctionTemplateInterpolationHistogram<
          TFCS1DFunction_HistogramInt32BinEdges,
          TFCS1DFunction_HistogramInt16BinEdges,
          uint32_t,
          float>
{
public:
  TFCS2DFunctionInt32Int16Int32InterpolationHistogram(TH2* h = nullptr)
      : TFCS2DFunctionTemplateInterpolationHistogram<
            TFCS1DFunction_HistogramInt32BinEdges,
            TFCS1DFunction_HistogramInt16BinEdges,
            uint32_t,
            float>(h) {};

  ClassDef(TFCS2DFunctionInt32Int16Int32InterpolationHistogram,
           1)  // TFCS2DFunctionInt32Int16Int32InterpolationHistogram
};

class TFCS2DFunctionInt32Int32Int8InterpolationHistogram
    : public TFCS2DFunctionTemplateInterpolationHistogram<
          TFCS1DFunction_HistogramInt32BinEdges,
          TFCS1DFunction_HistogramInt32BinEdges,
          uint8_t,
          float>
{
public:
  TFCS2DFunctionInt32Int32Int8InterpolationHistogram(TH2* h = nullptr)
      : TFCS2DFunctionTemplateInterpolationHistogram<
            TFCS1DFunction_HistogramInt32BinEdges,
            TFCS1DFunction_HistogramInt32BinEdges,
            uint8_t,
            float>(h) {};

  ClassDef(TFCS2DFunctionInt32Int32Int8InterpolationHistogram,
           1)  // TFCS2DFunctionInt32Int32Int8InterpolationHistogram
};
class TFCS2DFunctionInt32Int32Int16InterpolationHistogram
    : public TFCS2DFunctionTemplateInterpolationHistogram<
          TFCS1DFunction_HistogramInt32BinEdges,
          TFCS1DFunction_HistogramInt32BinEdges,
          uint16_t,
          float>
{
public:
  TFCS2DFunctionInt32Int32Int16InterpolationHistogram(TH2* h = nullptr)
      : TFCS2DFunctionTemplateInterpolationHistogram<
            TFCS1DFunction_HistogramInt32BinEdges,
            TFCS1DFunction_HistogramInt32BinEdges,
            uint16_t,
            float>(h) {};

  ClassDef(TFCS2DFunctionInt32Int32Int16InterpolationHistogram,
           1)  // TFCS2DFunctionInt32Int32Int16InterpolationHistogram
};
class TFCS2DFunctionInt32Int32Int32InterpolationHistogram
    : public TFCS2DFunctionTemplateInterpolationHistogram<
          TFCS1DFunction_HistogramInt32BinEdges,
          TFCS1DFunction_HistogramInt32BinEdges,
          uint32_t,
          float>
{
public:
  TFCS2DFunctionInt32Int32Int32InterpolationHistogram(TH2* h = nullptr)
      : TFCS2DFunctionTemplateInterpolationHistogram<
            TFCS1DFunction_HistogramInt32BinEdges,
            TFCS1DFunction_HistogramInt32BinEdges,
            uint32_t,
            float>(h) {};

  ClassDef(TFCS2DFunctionInt32Int32Int32InterpolationHistogram,
           1)  // TFCS2DFunctionInt32Int32Int32InterpolationHistogram
};

}  // namespace FastCaloSim::Core

// clang-format off
#if defined(__ROOTCLING__) && defined(__FastCaloSimStandAlone__)
#  pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateInterpolationHistogram \
          < FastCaloSim::Core::TFCS1DFunction_HistogramInt8BinEdges, \
      FastCaloSim::Core::TFCS1DFunction_HistogramInt8BinEdges, uint8_t, float> + ;
#  pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateInterpolationHistogram \
          < FastCaloSim::Core::TFCS1DFunction_HistogramInt8BinEdges, \
      FastCaloSim::Core::TFCS1DFunction_HistogramInt8BinEdges, uint16_t, float> + ;
#  pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateInterpolationHistogram \
          < FastCaloSim::Core::TFCS1DFunction_HistogramInt8BinEdges, \
      FastCaloSim::Core::TFCS1DFunction_HistogramInt8BinEdges, uint32_t, float> + ;
#  pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateInterpolationHistogram \
          < FastCaloSim::Core::TFCS1DFunction_HistogramInt8BinEdges, \
      FastCaloSim::Core::TFCS1DFunction_HistogramInt16BinEdges, uint8_t, float> + ;
#  pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateInterpolationHistogram \
          < FastCaloSim::Core::TFCS1DFunction_HistogramInt8BinEdges, \
      FastCaloSim::Core::TFCS1DFunction_HistogramInt16BinEdges, uint16_t, float> + ;
#  pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateInterpolationHistogram \
          < FastCaloSim::Core::TFCS1DFunction_HistogramInt8BinEdges, \
      FastCaloSim::Core::TFCS1DFunction_HistogramInt16BinEdges, uint32_t, float> + ;
#  pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateInterpolationHistogram \
          < FastCaloSim::Core::TFCS1DFunction_HistogramInt8BinEdges, \
      FastCaloSim::Core::TFCS1DFunction_HistogramInt32BinEdges, uint8_t, float> + ;
#  pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateInterpolationHistogram \
          < FastCaloSim::Core::TFCS1DFunction_HistogramInt8BinEdges, \
      FastCaloSim::Core::TFCS1DFunction_HistogramInt32BinEdges, uint16_t, float> + ;
#  pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateInterpolationHistogram \
          < FastCaloSim::Core::TFCS1DFunction_HistogramInt8BinEdges, \
      FastCaloSim::Core::TFCS1DFunction_HistogramInt32BinEdges, uint32_t, float> + ;
#  pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateInterpolationHistogram \
          < FastCaloSim::Core::TFCS1DFunction_HistogramInt16BinEdges, \
      FastCaloSim::Core::TFCS1DFunction_HistogramInt8BinEdges, uint8_t, float> + ;
#  pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateInterpolationHistogram \
          < FastCaloSim::Core::TFCS1DFunction_HistogramInt16BinEdges, \
      FastCaloSim::Core::TFCS1DFunction_HistogramInt8BinEdges, uint16_t, float> + ;
#  pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateInterpolationHistogram \
          < FastCaloSim::Core::TFCS1DFunction_HistogramInt16BinEdges, \
      FastCaloSim::Core::TFCS1DFunction_HistogramInt8BinEdges, uint32_t, float> + ;
#  pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateInterpolationHistogram \
          < FastCaloSim::Core::TFCS1DFunction_HistogramInt16BinEdges, \
      FastCaloSim::Core::TFCS1DFunction_HistogramInt16BinEdges, uint8_t, float> + ;
#  pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateInterpolationHistogram \
          < FastCaloSim::Core::TFCS1DFunction_HistogramInt16BinEdges, \
      FastCaloSim::Core::TFCS1DFunction_HistogramInt16BinEdges, uint16_t, float> + ;
#  pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateInterpolationHistogram \
          < FastCaloSim::Core::TFCS1DFunction_HistogramInt16BinEdges, \
      FastCaloSim::Core::TFCS1DFunction_HistogramInt16BinEdges, uint32_t, float> + ;
#  pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateInterpolationHistogram \
          < FastCaloSim::Core::TFCS1DFunction_HistogramInt16BinEdges, \
      FastCaloSim::Core::TFCS1DFunction_HistogramInt32BinEdges, uint8_t, float> + ;
#  pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateInterpolationHistogram \
          < FastCaloSim::Core::TFCS1DFunction_HistogramInt16BinEdges, \
      FastCaloSim::Core::TFCS1DFunction_HistogramInt32BinEdges, uint16_t, float> + ;
#  pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateInterpolationHistogram \
          < FastCaloSim::Core::TFCS1DFunction_HistogramInt16BinEdges, \
      FastCaloSim::Core::TFCS1DFunction_HistogramInt32BinEdges, uint32_t, float> + ;
#  pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateInterpolationHistogram \
          < FastCaloSim::Core::TFCS1DFunction_HistogramInt32BinEdges, \
      FastCaloSim::Core::TFCS1DFunction_HistogramInt8BinEdges, uint8_t, float> + ;
#  pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateInterpolationHistogram \
          < FastCaloSim::Core::TFCS1DFunction_HistogramInt32BinEdges, \
      FastCaloSim::Core::TFCS1DFunction_HistogramInt8BinEdges, uint16_t, float> + ;
#  pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateInterpolationHistogram \
          < FastCaloSim::Core::TFCS1DFunction_HistogramInt32BinEdges, \
      FastCaloSim::Core::TFCS1DFunction_HistogramInt8BinEdges, uint32_t, float> + ;
#  pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateInterpolationHistogram \
          < FastCaloSim::Core::TFCS1DFunction_HistogramInt32BinEdges, \
      FastCaloSim::Core::TFCS1DFunction_HistogramInt16BinEdges, uint8_t, float> + ;
#  pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateInterpolationHistogram \
          < FastCaloSim::Core::TFCS1DFunction_HistogramInt32BinEdges, \
      FastCaloSim::Core::TFCS1DFunction_HistogramInt16BinEdges, uint16_t, float> + ;
#  pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateInterpolationHistogram \
          < FastCaloSim::Core::TFCS1DFunction_HistogramInt32BinEdges, \
      FastCaloSim::Core::TFCS1DFunction_HistogramInt16BinEdges, uint32_t, float> + ;
#  pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateInterpolationHistogram \
          < FastCaloSim::Core::TFCS1DFunction_HistogramInt32BinEdges, \
      FastCaloSim::Core::TFCS1DFunction_HistogramInt32BinEdges, uint8_t, float> + ;
#  pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateInterpolationHistogram \
          < FastCaloSim::Core::TFCS1DFunction_HistogramInt32BinEdges, \
      FastCaloSim::Core::TFCS1DFunction_HistogramInt32BinEdges, uint16_t, float> + ;
#  pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateInterpolationHistogram \
          < FastCaloSim::Core::TFCS1DFunction_HistogramInt32BinEdges, \
      FastCaloSim::Core::TFCS1DFunction_HistogramInt32BinEdges, uint32_t, float> + ;

#  pragma link C++ class FastCaloSim::Core::TFCS2DFunctionInt8Int8Int8InterpolationHistogram + ;
#  pragma link C++ class FastCaloSim::Core::TFCS2DFunctionInt8Int8Int16InterpolationHistogram + ;
#  pragma link C++ class FastCaloSim::Core::TFCS2DFunctionInt8Int8Int32InterpolationHistogram + ;
#  pragma link C++ class FastCaloSim::Core::TFCS2DFunctionInt8Int16Int8InterpolationHistogram + ;
#  pragma link C++ class FastCaloSim::Core::TFCS2DFunctionInt8Int16Int16InterpolationHistogram + ;
#  pragma link C++ class FastCaloSim::Core::TFCS2DFunctionInt8Int16Int32InterpolationHistogram + ;
#  pragma link C++ class FastCaloSim::Core::TFCS2DFunctionInt8Int32Int8InterpolationHistogram + ;
#  pragma link C++ class FastCaloSim::Core::TFCS2DFunctionInt8Int32Int16InterpolationHistogram + ;
#  pragma link C++ class FastCaloSim::Core::TFCS2DFunctionInt8Int32Int32InterpolationHistogram + ;
#  pragma link C++ class FastCaloSim::Core::TFCS2DFunctionInt16Int8Int8InterpolationHistogram + ;
#  pragma link C++ class FastCaloSim::Core::TFCS2DFunctionInt16Int8Int16InterpolationHistogram + ;
#  pragma link C++ class FastCaloSim::Core::TFCS2DFunctionInt16Int8Int32InterpolationHistogram + ;
#  pragma link C++ class FastCaloSim::Core::TFCS2DFunctionInt16Int16Int8InterpolationHistogram + ;
#  pragma link C++ class FastCaloSim::Core::TFCS2DFunctionInt16Int16Int16InterpolationHistogram + ;
#  pragma link C++ class FastCaloSim::Core::TFCS2DFunctionInt16Int16Int32InterpolationHistogram + ;
#  pragma link C++ class FastCaloSim::Core::TFCS2DFunctionInt16Int32Int8InterpolationHistogram + ;
#  pragma link C++ class FastCaloSim::Core::TFCS2DFunctionInt16Int32Int16InterpolationHistogram + ;
#  pragma link C++ class FastCaloSim::Core::TFCS2DFunctionInt16Int32Int32InterpolationHistogram + ;
#  pragma link C++ class FastCaloSim::Core::TFCS2DFunctionInt32Int8Int8InterpolationHistogram + ;
#  pragma link C++ class FastCaloSim::Core::TFCS2DFunctionInt32Int8Int16InterpolationHistogram + ;
#  pragma link C++ class FastCaloSim::Core::TFCS2DFunctionInt32Int8Int32InterpolationHistogram + ;
#  pragma link C++ class FastCaloSim::Core::TFCS2DFunctionInt32Int16Int8InterpolationHistogram + ;
#  pragma link C++ class FastCaloSim::Core::TFCS2DFunctionInt32Int16Int16InterpolationHistogram + ;
#  pragma link C++ class FastCaloSim::Core::TFCS2DFunctionInt32Int16Int32InterpolationHistogram + ;
#  pragma link C++ class FastCaloSim::Core::TFCS2DFunctionInt32Int32Int8InterpolationHistogram + ;
#  pragma link C++ class FastCaloSim::Core::TFCS2DFunctionInt32Int32Int16InterpolationHistogram + ;
#  pragma link C++ class FastCaloSim::Core::TFCS2DFunctionInt32Int32Int32InterpolationHistogram + ;

#endif
// clang-format on

#endif
