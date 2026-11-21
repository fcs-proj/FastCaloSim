// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#pragma once

#include "FastCaloSim/Core/TFCS1DFunctionTemplateHistogram.h"

template<typename Txvec, typename Ty, typename Trandom = float>
class TFCS1DFunctionTemplateInterpolationHistogram
    : public TFCS1DFunctionTemplateHistogram<Txvec, Ty, Trandom>
{
public:
  using TFCS1DFunctionTemplateHistogram<Txvec, Ty, Trandom>::
      TFCS1DFunctionTemplateHistogram;

  using TFCS1DFunction::rnd_to_fct;
  using TFCS1DFunctionTemplateHistogram<Txvec, Ty, Trandom>::get_nbins;
  using TFCS1DFunctionTemplateHistogram<Txvec, Ty, Trandom>::m_HistoContents;
  using TFCS1DFunctionTemplateHistogram<Txvec, Ty, Trandom>::m_HistoBorders;

  /// Function gets random number rnd in the range [0,1) as argument
  /// and returns function value according to a histogram distribution.
  /// A linear interpolation is done within each bin
  auto rnd_to_fct(double rnd) const -> double override
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

    // Slope of the density inside the selected bin. The density in a bin is
    // dfrac/l. Its change to the next and to the previous bin, divided by the
    // distance between the bin centers (l+lnext)/2 and (l+lprev)/2, is averaged
    // and multiplied by l to get the change across this bin. Normalized to the
    // density in this bin, dfrac/l, it is the slope m in units of the bin width
    // that position_lin expects.
    Trandom m;
    if (dfrac > 0) {
      Trandom mnext = dfracnext / lnext - dfrac / l;
      mnext /= l + lnext;
      Trandom mprev = dfrac / l - dfracprev / lprev;
      mprev /= l + lprev;
      m = (mnext + mprev) * l * l / dfrac;
    } else
      m = 0;

    return m_HistoBorders.position_lin(ibin, m, residual_rnd);
  }

  ClassDefOverride(TFCS1DFunctionTemplateInterpolationHistogram,
                   1)  // TFCS1DFunctionTemplateInterpolationHistogram
};

class TFCS1DFunctionInt8Int8InterpolationHistogram
    : public TFCS1DFunctionTemplateInterpolationHistogram<
          TFCS1DFunction_HistogramInt8BinEdges,
          uint8_t>
{
public:
  using TFCS1DFunctionTemplateInterpolationHistogram::
      TFCS1DFunctionTemplateInterpolationHistogram;

  ClassDefOverride(TFCS1DFunctionInt8Int8InterpolationHistogram,
                   1)  // TFCS1DFunctionInt8Int8InterpolationHistogram
};

class TFCS1DFunctionInt8Int16InterpolationHistogram
    : public TFCS1DFunctionTemplateInterpolationHistogram<
          TFCS1DFunction_HistogramInt8BinEdges,
          uint16_t>
{
public:
  using TFCS1DFunctionTemplateInterpolationHistogram::
      TFCS1DFunctionTemplateInterpolationHistogram;

  ClassDefOverride(TFCS1DFunctionInt8Int16InterpolationHistogram,
                   1)  // TFCS1DFunctionInt8Int16InterpolationHistogram
};
