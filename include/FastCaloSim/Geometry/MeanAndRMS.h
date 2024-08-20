/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

#ifndef MeanAndRMS_h
#define MeanAndRMS_h

#include <math.h>

class MeanAndRMS
{
public:
  MeanAndRMS()
      : m_w(0)
      , m_wx(0)
      , m_wx2(0) {};
  MeanAndRMS(const double xadd, const double weight = 1)
      : m_w(weight)
      , m_wx(weight * xadd)
      , m_wx2(weight * xadd * xadd) {};
  MeanAndRMS(const MeanAndRMS& ref)
      : m_w(ref.m_w)
      , m_wx(ref.m_wx)
      , m_wx2(ref.m_wx2) {};

  MeanAndRMS& add(double xadd, double weight = 1)
  {
    m_wx += weight * xadd;
    m_wx2 += weight * xadd * xadd;
    m_w += weight;
    return *this;
  };
  MeanAndRMS& operator+=(double xadd) { return add(xadd); };
  MeanAndRMS& operator-=(double xadd) { return add(-xadd); };
  MeanAndRMS& operator=(const MeanAndRMS& ref)
  {
    m_w = ref.m_w;
    m_wx = ref.m_wx;
    m_wx2 = ref.m_wx2;
    return *this;
  };

  double sum_weight() const { return m_w; };
  double mean() const
  {
    if (m_w != 0)
      return m_wx / m_w;
    else
      return 0;
  };
  double mean2() const
  {
    double x = mean();
    return x * x;
  };
  double rms2() const
  {
    if (m_w != 0)
      return m_wx2 / m_w - mean2();
    else
      return 0;
  };
  double rms() const
  {
    double r2 = rms2();
    if (r2 >= 0)
      return sqrt(r2);
    else
      return 0;
  };
  double mean_error() const
  {
    if (m_w > 0)
      return rms() / sqrt(m_w);
    else
      return 0;
  };
  double rms_error() const
  {
    if (m_w > 0)
      return rms() / sqrt(2 * m_w);
    else
      return 0;
  };

  operator double() const { return mean(); }

protected:
  double m_w;
  double m_wx, m_wx2;
};

#endif
