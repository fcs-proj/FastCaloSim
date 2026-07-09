// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#ifndef ISF_FASTCALOSIMEVENT_TFCS1DFunctionSpline_h
#define ISF_FASTCALOSIMEVENT_TFCS1DFunctionSpline_h

#include <vector>

#include <TSpline.h>

#include "FastCaloSim/Core/TFCS1DFunction.h"

class TH1;
class TFCS1DFunctionInt32Histogram;

class TFCS1DFunctionSpline : public TFCS1DFunction
//    ^ heavy use of gRandom, but class is not used anywhere
{
public:
  TFCS1DFunctionSpline(TH1* hist = nullptr,
                       double maxdevgoal = 0.01,
                       double maxeffsiggoal = 3,
                       int maxnp = 20)
  {
    if (hist)
      Initialize(hist, maxdevgoal, maxeffsiggoal, maxnp);
  };
  ~TFCS1DFunctionSpline() override = default;

  static auto get_maxdev(const TH1* hist,
                         const TSpline3& sp,
                         double& maxeffsig,
                         double& p_maxdev,
                         double& p_maxeffsig,
                         int ntoy = 10000) -> double;

  auto Initialize(TH1* hist,
                  double maxdevgoal = 0.01,
                  double maxeffsiggoal = 3,
                  int maxnp = 20) -> double;

  auto InitializeFromSpline(TH1* hist,
                            const TSpline3& sp,
                            double maxdevgoal = 0.01,
                            double maxeffsiggoal = 3) -> double;
  auto InitializeEqualDistance(TH1* hist,
                               double maxdevgoal = 0.01,
                               double maxeffsiggoal = 3,
                               int nsplinepoints = 5) -> double;
  auto InitializeEqualProbability(TH1* hist,
                                  double maxdevgoal = 0.01,
                                  double maxeffsiggoal = 3,
                                  int nsplinepoints = 5) -> double;

  using TFCS1DFunction::rnd_to_fct;

  /// Function gets random number rnd in the range [0,1) as argument
  /// and returns function value according to a histogram distribution
  auto rnd_to_fct(double rnd) const -> double override;

  auto spline() const -> const TSpline3& { return m_spline; };
  auto spline() -> TSpline3& { return m_spline; };

protected:
  static auto optimize(TSpline3& sp_best,
                       std::vector<double>& nprop,
                       const TH1* hist,
                       TFCS1DFunctionInt32Histogram& hist_fct,
                       double maxdevgoal = 0.01,
                       double maxeffsiggoal = 3) -> double;

  TSpline3 m_spline;

private:
  ClassDef(TFCS1DFunctionSpline, 1)  // TFCS1DFunctionSpline
};

#endif
