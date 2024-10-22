// Copyright (c) 2024 CERN for the benefit of the FastCaloSim project

#ifndef ISF_FASTCALOSIMEVENT_TFCS1DFunctionTemplateInterpolationHistogram_h
#define ISF_FASTCALOSIMEVENT_TFCS1DFunctionTemplateInterpolationHistogram_h

#include "FastCaloSim/Core/TFCS1DFunctionTemplateHistogram.h"

template<typename Txvec, typename Ty, typename Trandom = float>
class TFCS1DFunctionTemplateInterpolationHistogram
    : public TFCS1DFunctionTemplateHistogram<Txvec, Ty, Trandom>
{
public:
  TFCS1DFunctionTemplateInterpolationHistogram(TH1* hist = nullptr)
      : TFCS1DFunctionTemplateHistogram<Txvec, Ty, Trandom>(hist) {};

  using TFCS1DFunction::rnd_to_fct;
  using TFCS1DFunctionTemplateHistogram<Txvec, Ty, Trandom>::get_nbins;
  using TFCS1DFunctionTemplateHistogram<Txvec, Ty, Trandom>::m_HistoContents;
  using TFCS1DFunctionTemplateHistogram<Txvec, Ty, Trandom>::m_HistoBorders;

  /// Function gets random number rnd in the range [0,1) as argument
  /// and returns function value according to a histogram distribution.
  /// A linear interpolation is done within each bin
  virtual double rnd_to_fct(double rnd) const
  {
    size_t nbins = get_nbins();
    if (nbins == 0)
      return 0;
    Trandom residual_rnd;
    size_t ibin = m_HistoContents.get_bin(rnd, residual_rnd);

    Trandom frac = m_HistoContents.get_fraction(ibin);

    Trandom dfrac;
    auto l = m_HistoBorders.GetBinLength(ibin);
    Trandom dfracprev;
    auto lprev = l;
    if (ibin > 0) {
      Trandom fracprev = m_HistoContents.get_fraction(ibin - 1);
      dfrac = frac - fracprev;
      if (ibin > 1)
        dfracprev = fracprev - m_HistoContents.get_fraction(ibin - 2);
      else
        dfracprev = fracprev;
      lprev = m_HistoBorders.GetBinLength(ibin - 1);
    } else {
      dfrac = frac;
      dfracprev = dfrac;
    }

    Trandom dfracnext;
    auto lnext = l;
    if (ibin < nbins - 1) {
      dfracnext = m_HistoContents.get_fraction(ibin + 1) - frac;
      lnext = m_HistoBorders.GetBinLength(ibin + 1);
    } else {
      dfracnext = dfrac;
    }

    Trandom m;
    if (dfrac > 0) {
      /*
      //approximately readable version of the code below
      //density in bin is dfracX/lX
      Trandom mnext=dfracnext/lnext-dfrac/l; //change in density to next bin
      mnext/=0.5*(l+lnext); //normalize change in density to distance between
      bins Trandom mprev=dfrac/l-dfracprev/lprev; //change in density to
      previous bin mprev/=0.5*(l+lprev); //normalize change in density to
      distance between bins m=0.5*(mnext+mprev); //average the two m*=l; //m is
      needed inside this bin, so multiply by size of this bin (matches
      normalization to distance between bins above) m/=dfrac/l; //finally, m is
      the slope in this bin, given that this bin was already selected. So
      overall normalize to density inside this bin
      */
      Trandom mnext = dfracnext / lnext - dfrac / l;
      mnext /= l + lnext;
      Trandom mprev = dfrac / l - dfracprev / lprev;
      mprev /= l + lprev;
      m = (mnext + mprev) * l * l / dfrac;
    } else
      m = 0;

    // ATH_MSG_INFO( fx="<<m_HistoBorders.GetBinLowEdge(ibin)<<"
    // frac="<<m_HistoContents.get_fraction(ibin)<<" dfracprev="<<dfracprev<<"
    // dfrac="<<dfrac<<" dfracnext="<<dfracnext<<"
    // dfracprev-dfrac="<<dfracprev-dfrac<<"
    // dfracnext-dfrac="<<dfracnext-dfrac<<" m="<<m<<"
    // residual_rnd="<<residual_rnd);
    return m_HistoBorders.position_lin(ibin, m, residual_rnd);
  }

  ClassDef(TFCS1DFunctionTemplateInterpolationHistogram,
           1)  // TFCS1DFunctionTemplateInterpolationHistogram
};

class TFCS1DFunctionInt8Int8InterpolationHistogram
    : public TFCS1DFunctionTemplateInterpolationHistogram<
          TFCS1DFunction_HistogramInt8BinEdges,
          uint8_t,
          float>
{
public:
  TFCS1DFunctionInt8Int8InterpolationHistogram(TH1* h = nullptr)
      : TFCS1DFunctionTemplateInterpolationHistogram<
          TFCS1DFunction_HistogramInt8BinEdges,
          uint8_t,
          float>(h) {};

  ClassDef(TFCS1DFunctionInt8Int8InterpolationHistogram,
           1)  // TFCS1DFunctionInt8Int8InterpolationHistogram
};

class TFCS1DFunctionInt8Int16InterpolationHistogram
    : public TFCS1DFunctionTemplateInterpolationHistogram<
          TFCS1DFunction_HistogramInt8BinEdges,
          uint16_t,
          float>
{
public:
  TFCS1DFunctionInt8Int16InterpolationHistogram(TH1* h = nullptr)
      : TFCS1DFunctionTemplateInterpolationHistogram<
          TFCS1DFunction_HistogramInt8BinEdges,
          uint16_t,
          float>(h) {};

  ClassDef(TFCS1DFunctionInt8Int16InterpolationHistogram,
           1)  // TFCS1DFunctionInt8Int16InterpolationHistogram
};

class TFCS1DFunctionInt16Int16InterpolationHistogram
    : public TFCS1DFunctionTemplateInterpolationHistogram<
          TFCS1DFunction_HistogramInt16BinEdges,
          uint16_t,
          float>
{
public:
  TFCS1DFunctionInt16Int16InterpolationHistogram(TH1* h = nullptr)
      : TFCS1DFunctionTemplateInterpolationHistogram<
          TFCS1DFunction_HistogramInt16BinEdges,
          uint16_t,
          float>(h) {};

  ClassDef(TFCS1DFunctionInt16Int16InterpolationHistogram,
           1)  // TFCS1DFunctionInt16Int16InterpolationHistogram
};

class TFCS1DFunctionInt16Int32InterpolationHistogram
    : public TFCS1DFunctionTemplateInterpolationHistogram<
          TFCS1DFunction_HistogramInt16BinEdges,
          uint32_t,
          float>
{
public:
  TFCS1DFunctionInt16Int32InterpolationHistogram(TH1* h = nullptr)
      : TFCS1DFunctionTemplateInterpolationHistogram<
          TFCS1DFunction_HistogramInt16BinEdges,
          uint32_t,
          float>(h) {};

  ClassDef(TFCS1DFunctionInt16Int32InterpolationHistogram,
           1)  // TFCS1DFunctionInt16Int32InterpolationHistogram
};

// clang-format off
#if defined(__ROOTCLING__) && defined(__FastCaloSimStandAlone__)
#  pragma link C++ class TFCS1DFunctionTemplateInterpolationHistogram \
          < TFCS1DFunction_HistogramInt8BinEdges, \
      uint8_t, float> + ;
#  pragma link C++ class TFCS1DFunctionTemplateInterpolationHistogram \
          < TFCS1DFunction_HistogramInt8BinEdges, \
      uint16_t, float> + ;
#  pragma link C++ class TFCS1DFunctionTemplateInterpolationHistogram \
          < TFCS1DFunction_HistogramInt16BinEdges, \
      uint16_t, float> + ;
#  pragma link C++ class TFCS1DFunctionTemplateInterpolationHistogram \
          < TFCS1DFunction_HistogramInt16BinEdges, \
      uint32_t, float> + ;

#  pragma link C++ class TFCS1DFunctionInt8Int8InterpolationHistogram + ;
#  pragma link C++ class TFCS1DFunctionInt8Int16InterpolationHistogram + ;
#  pragma link C++ class TFCS1DFunctionInt16Int16InterpolationHistogram + ;
#  pragma link C++ class TFCS1DFunctionInt16Int32InterpolationHistogram + ;

#endif
// clang-format on

#endif
