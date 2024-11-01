// Copyright (c) 2024 CERN for the benefit of the FastCaloSim project

#include <iostream>

#include "FastCaloSim/Core/TFCS1DFunctionSpline.h"

#include "FastCaloSim/Core/TFCS1DFunctionInt32Histogram.h"
#include "TH1.h"
#include "TMath.h"
#include "TRandom.h"

//=============================================
//======= TFCS1DFunctionSpline =========
//=============================================

double TFCS1DFunctionSpline::Initialize(TH1* hist,
                                        double maxdevgoal,
                                        double maxeffsiggoal,
                                        int maxnp)
{
  double max_penalty_best = -1;
  TSpline3 sp_best;
  for (int np = 3; np <= maxnp; ++np) {
    FCS_MSG_INFO("========== Spline #=" << np << " ==============");
    double max_penalty;
    if (max_penalty_best > 0) {
      max_penalty =
          InitializeFromSpline(hist, sp_best, maxdevgoal, maxeffsiggoal);
      if (max_penalty_best < 0 || max_penalty < max_penalty_best) {
        max_penalty_best = max_penalty;
        sp_best = m_spline;
      }
    }

    max_penalty = InitializeEqualDistance(hist, maxdevgoal, maxeffsiggoal, np);
    if (max_penalty_best < 0 || max_penalty < max_penalty_best) {
      max_penalty_best = max_penalty;
      sp_best = m_spline;
    }

    max_penalty =
        InitializeEqualProbability(hist, maxdevgoal, maxeffsiggoal, np);
    if (max_penalty_best < 0 || max_penalty < max_penalty_best) {
      max_penalty_best = max_penalty;
      sp_best = m_spline;
    }

    FCS_MSG_INFO("========== Spline #=" << np << " max_penalty_best="
                                        << max_penalty_best
                                        << " ==============");
    FCS_MSG(INFO) << "==== Best spline init | ";
    for (int i = 0; i < sp_best.GetNp(); ++i) {
      double p, x;
      sp_best.GetKnot(i, p, x);
      FCS_MSG(INFO) << i << " : p=" << p << " x=" << x << " ; ";
    }
    FCS_MSG(INFO) << " =====" << END_FCS_MSG(INFO);

    if (max_penalty_best < 2)
      break;
  }
  m_spline = sp_best;

  return max_penalty_best;
}

double TFCS1DFunctionSpline::InitializeFromSpline(TH1* hist,
                                                  const TSpline3& sp,
                                                  double maxdevgoal,
                                                  double maxeffsiggoal)
{
  TFCS1DFunctionInt32Histogram hist_fct(hist);

  double maxeffsig;
  double p_maxdev;
  double p_maxeffsig;
  double maxdev = get_maxdev(hist, sp, maxeffsig, p_maxdev, p_maxeffsig);
  double p_improve;
  if (maxdev / maxdevgoal > maxeffsig / maxeffsiggoal)
    p_improve = p_maxdev;
  else
    p_improve = p_maxeffsig;

  int nsplinepoints = sp.GetNp();
  std::vector<double> nprop(nsplinepoints + 1);
  int ind = 0;
  FCS_MSG(INFO) << "Spline init p_improve=" << p_improve << " | ";
  for (int i = 0; i < nsplinepoints; ++i) {
    double p, x;
    sp.GetKnot(i, p, x);
    if (i == 0 && p_improve < p) {
      nprop[ind] = (0 + p) / 2;
      FCS_MSG(INFO) << ind << " : pi=" << nprop[ind] << " ; ";
      ++ind;
    }

    nprop[ind] = p;
    FCS_MSG(INFO) << ind << " : p=" << nprop[ind] << " ; ";
    ++ind;

    if (i == nsplinepoints - 1 && p_improve > p) {
      nprop[ind] = (1 + p) / 2;
      FCS_MSG(INFO) << ind << " : pi=" << nprop[ind] << " ; ";
      ++ind;
    }
    if (i < nsplinepoints - 1) {
      double pn, xn;
      sp.GetKnot(i + 1, pn, xn);
      if (p_improve > p && p_improve < pn) {
        nprop[ind] = (p + pn) / 2;
        FCS_MSG(INFO) << ind << " : pi=" << nprop[ind] << " ; ";
        ++ind;
      }
    }
  }
  FCS_MSG(INFO) << END_FCS_MSG(INFO);
  nsplinepoints = ind;
  nprop.resize(nsplinepoints);

  double max_penalty =
      optimize(m_spline, nprop, hist, hist_fct, maxdevgoal, maxeffsiggoal);
  maxdev = get_maxdev(hist, m_spline, maxeffsig, p_maxdev, p_maxeffsig);
  FCS_MSG_INFO("Spline init spline #=" << nsplinepoints << " : maxdev="
                                       << maxdev << " p_maxdev=" << p_maxdev
                                       << " maxeffsig=" << maxeffsig
                                       << " p_maxeffsig=" << p_maxeffsig
                                       << " max_penalty=" << max_penalty);
  FCS_MSG(INFO) << "  ";
  for (int i = 0; i < m_spline.GetNp(); ++i) {
    double p, x;
    m_spline.GetKnot(i, p, x);
    FCS_MSG(INFO) << i << " : p=" << p << " x=" << x << " ; ";
  }
  FCS_MSG(INFO) << END_FCS_MSG(INFO);
  return max_penalty;
}

double TFCS1DFunctionSpline::InitializeEqualDistance(TH1* hist,
                                                     double maxdevgoal,
                                                     double maxeffsiggoal,
                                                     int nsplinepoints)
{
  TFCS1DFunctionInt32Histogram hist_fct(hist);

  double xmin = 0;
  double xmax = 0;
  for (int i = 1; i <= hist->GetNbinsX(); i++) {
    xmin = hist->GetXaxis()->GetBinLowEdge(i);
    if (hist->GetBinContent(i) > 0)
      break;
  }
  for (int i = hist->GetNbinsX(); i >= 1; i--) {
    xmax = hist->GetXaxis()->GetBinUpEdge(i);
    if (hist->GetBinContent(i) > 0)
      break;
  }
  // FCS_MSG_INFO("xmin="<<xmin<<" xmax="<<xmax);

  double dx = (xmax - xmin) / (nsplinepoints - 1);

  std::vector<double> nprop(nsplinepoints);
  std::vector<double> nx(nsplinepoints);
  nprop[0] = 0;
  nx[0] = hist_fct.rnd_to_fct(nprop[0]);
  // FCS_MSG_INFO(0<<" p="<<nprop[0]<<" x="<<nx[0]);
  for (int i = 1; i < nsplinepoints; ++i) {
    nx[i] = xmin + i * dx;
    double p_min = 0;
    double p_max = 1;
    double p_test;
    double tx;
    do {
      p_test = 0.5 * (p_min + p_max);
      tx = hist_fct.rnd_to_fct(p_test);
      if (nx[i] < tx)
        p_max = p_test;
      else
        p_min = p_test;
      if ((p_max - p_min) < 0.0000001)
        break;
    } while (TMath::Abs(tx - nx[i]) > dx / 10);
    // FCS_MSG_INFO(i<<" p="<<p_test<<" x="<<tx);
    nprop[i] = p_test;
  }

  double max_penalty =
      optimize(m_spline, nprop, hist, hist_fct, maxdevgoal, maxeffsiggoal);
  double maxeffsig;
  double p_maxdev;
  double p_maxeffsig;
  double maxdev = get_maxdev(hist, m_spline, maxeffsig, p_maxdev, p_maxeffsig);
  FCS_MSG_INFO("Spline init equ. dist. #=" << nsplinepoints << " : maxdev="
                                           << maxdev << " p_maxdev=" << p_maxdev
                                           << " maxeffsig=" << maxeffsig
                                           << " p_maxeffsig=" << p_maxeffsig
                                           << " max_penalty=" << max_penalty);
  FCS_MSG(INFO) << "  ";
  for (int i = 0; i < m_spline.GetNp(); ++i) {
    double p, x;
    m_spline.GetKnot(i, p, x);
    FCS_MSG(INFO) << i << " : p=" << p << " x=" << x << " ; ";
  }
  FCS_MSG(INFO) << END_FCS_MSG(INFO);
  return max_penalty;
}

double TFCS1DFunctionSpline::InitializeEqualProbability(TH1* hist,
                                                        double maxdevgoal,
                                                        double maxeffsiggoal,
                                                        int nsplinepoints)
{
  TFCS1DFunctionInt32Histogram hist_fct(hist);

  double dprop = 1.0 / (nsplinepoints - 1);
  std::vector<double> nprop(nsplinepoints);
  for (int i = 0; i < nsplinepoints; ++i) {
    nprop[i] = i * dprop;
  }

  double max_penalty =
      optimize(m_spline, nprop, hist, hist_fct, maxdevgoal, maxeffsiggoal);
  double maxeffsig;
  double p_maxdev;
  double p_maxeffsig;
  double maxdev = get_maxdev(hist, m_spline, maxeffsig, p_maxdev, p_maxeffsig);
  FCS_MSG_INFO("Spline init equ. prob. #=" << nsplinepoints << " : maxdev="
                                           << maxdev << " p_maxdev=" << p_maxdev
                                           << " maxeffsig=" << maxeffsig
                                           << " p_maxeffsig=" << p_maxeffsig
                                           << " max_penalty=" << max_penalty);
  FCS_MSG(INFO) << "  ";
  for (int i = 0; i < m_spline.GetNp(); ++i) {
    double p, x;
    m_spline.GetKnot(i, p, x);
    FCS_MSG(INFO) << i << " : p=" << p << " x=" << x << " ; ";
  }
  FCS_MSG(INFO) << END_FCS_MSG(INFO);
  return max_penalty;
}

double TFCS1DFunctionSpline::optimize(TSpline3& sp_best,
                                      std::vector<double>& nprop,
                                      const TH1* hist,
                                      TFCS1DFunctionInt32Histogram& hist_fct,
                                      double maxdevgoal,
                                      double maxeffsiggoal)
{
  int nsplinepoints = (int)nprop.size();
  // double xmin=hist->GetXaxis()->GetXmin();
  // double xmax=hist->GetXaxis()->GetXmax();
  std::vector<double> nx(nsplinepoints);
  std::vector<double> nprop_try = nprop;
  double max_penalty = -1;
  double p_gotobest = 0.5 / nsplinepoints;
  int ntry = 200 / p_gotobest;
  // int itrytot=0;
  for (double dproploop = 0.4 / nsplinepoints; dproploop > 0.02 / nsplinepoints;
       dproploop /= 2)
  {
    double dprop = dproploop;
    int n_nogain = 0;
    for (int itry = 0; itry < ntry; ++itry) {
      // itrytot++;
      for (int i = 0; i < nsplinepoints; ++i) {
        nx[i] = hist_fct.rnd_to_fct(nprop_try[i]);
      }
      TSpline3 sp(
          "1Dspline", nprop_try.data(), nx.data(), nsplinepoints, "b2e2", 0, 0);
      double maxeffsig;
      double p_maxdev;
      double p_maxeffsig;
      double maxdev = get_maxdev(hist, sp, maxeffsig, p_maxdev, p_maxeffsig);
      double penalty = maxdev / maxdevgoal + maxeffsig / maxeffsiggoal;
      if (max_penalty < 0 || penalty < max_penalty) {
        max_penalty = penalty;
        nprop = nprop_try;
        sp_best = sp;
        /*
        FCS_MSG(INFO) <<"#="<<nsplinepoints<<" try="<<itrytot-1<<" | ";
        for(int i=0;i<nsplinepoints;++i) {
          FCS_MSG(INFO) <<i<<":p="<<nprop_try[i]<<" x="<<nx[i]<<" ; ";
        }
        FCS_MSG(INFO) <<"new maxdev="<<maxdev<<" maxeffsig="<<maxeffsig<<"
        max_penalty="<<max_penalty<<END_FCS_MSG(INFO);
        */
        n_nogain = 0;
      } else {
        if (gRandom->Rndm() < p_gotobest)
          nprop_try = nprop;
        ++n_nogain;
        // allow ~20 times retrying from the best found spline before aborting
        if (n_nogain > 20 / p_gotobest)
          break;
      }
      int ip = 1 + gRandom->Rndm() * (nsplinepoints - 1);
      if (ip > nsplinepoints - 1)
        ip = nsplinepoints - 1;
      double d = dprop;
      if (gRandom->Rndm() > 0.5)
        d = -dprop;
      double nprop_new = nprop_try[ip] + d;
      if (ip > 0)
        if (nprop_try[ip - 1] + dprop / 2 > nprop_new) {
          nprop_new = nprop_try[ip];
          dprop /= 2;
        }
      if (ip < nsplinepoints - 1)
        if (nprop_new > nprop_try[ip + 1] - dprop / 2) {
          nprop_new = nprop_try[ip];
          dprop /= 2;
        }
      if (nprop_new < 0) {
        nprop_new = 0;
        dprop /= 2;
      }
      if (nprop_new > 1) {
        nprop_new = 1;
        dprop /= 2;
      }
      nprop_try[ip] = nprop_new;
    }
    nprop_try = nprop;
  }
  return max_penalty;
}

double TFCS1DFunctionSpline::get_maxdev(const TH1* hist,
                                        const TSpline3& sp,
                                        double& maxeffsig,
                                        double& p_maxdev,
                                        double& p_maxeffsig,
                                        int ntoy)
{
  double maxdev = 0;
  maxeffsig = 0;

  TH1* hist_clone = (TH1*)hist->Clone("hist_clone");
  hist_clone->SetDirectory(nullptr);
  hist_clone->Reset();
  double interr = 0;
  double integral = hist->IntegralAndError(1, hist->GetNbinsX(), interr);
  double effN = integral / interr;
  effN *= effN;
  // FCS_MSG_INFO("integral="<<integral<<" +- "<<interr<<"
  // relerr="<<interr/integral); FCS_MSG_INFO("effN="<<effN<<" +-
  // "<<TMath::Sqrt(effN)<<" relerr="<<1/TMath::Sqrt(effN));
  double toyweight = 1.0 / ntoy;
  for (int itoy = 0; itoy < ntoy; ++itoy) {
    double prop = itoy * toyweight;
    hist_clone->Fill(sp.Eval(prop), toyweight);
  }

  double int1 = 0;
  double int2 = 0;
  for (int i = 0; i <= hist->GetNbinsX() + 1; i++) {
    int1 += hist->GetBinContent(i) / integral;
    int2 += hist_clone->GetBinContent(i);
    double val = TMath::Abs(int1 - int2);
    if (val > maxdev) {
      maxdev = val;
      p_maxdev = int1;
    }

    // now view the normalized integral as selection efficiency from a total
    // sample of sizze effN
    double int1err = TMath::Sqrt(int1 * (1 - int1) / effN);
    double valsig = 0;
    if (int1err > 0)
      valsig = val / int1err;
    if (valsig > maxeffsig) {
      maxeffsig = valsig;
      p_maxeffsig = int1;
    }

    // FCS_MSG_INFO(i<<": diff="<<int1-int2<<" sig(diff)="<<valsig<<"
    // int1="<<int1<<" +- "<<int1err<<" int2="<<int2<<" maxdev="<<maxdev<<"
    // maxeffsig="<<maxeffsig);
  }

  delete hist_clone;

  return maxdev;
}

double TFCS1DFunctionSpline::rnd_to_fct(double rnd) const
{
  return m_spline.Eval(rnd);
}
