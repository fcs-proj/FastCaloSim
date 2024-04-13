/*
  Copyright (C) 2002-2020 CERN for the benefit of the ATLAS collaboration
*/

#include "FastCaloSim/TFCSEnergyInterpolationSpline.h"
#include "FastCaloSim/TFCSSimulationState.h"
#include "FastCaloSim/TFCSTruthState.h"
#include "FastCaloSim/TFCSExtrapolationState.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TAxis.h"
#include <iostream>
#include <vector>

namespace Gaudi {
namespace Units {
constexpr double megaelectronvolt = 1.;
constexpr double kiloelectronvolt = 1.e-3 * megaelectronvolt;
constexpr double keV = kiloelectronvolt;
} // namespace Units
} // namespace Gaudi

//=============================================
//======= TFCSEnergyInterpolation =========
//=============================================

TFCSEnergyInterpolationSpline::TFCSEnergyInterpolationSpline(const char *name,
                                                             const char *title)
    : TFCSParametrization(name, title) {}

void TFCSEnergyInterpolationSpline::InitFromArrayInLogEkin(
    Int_t np, Double_t logEkin[], Double_t response[], const char *opt,
    Double_t valbeg, Double_t valend) {
  TSpline3 initspline(GetName(), logEkin, response, np, opt, valbeg, valend);
  m_spline = initspline;
}

void TFCSEnergyInterpolationSpline::InitFromArrayInEkin(
    Int_t np, Double_t Ekin[], Double_t response[], const char *opt,
    Double_t valbeg, Double_t valend) {
  std::vector<Double_t> logEkin(np);
  for (int i = 0; i < np; ++i)
    logEkin[i] = TMath::Log(Ekin[i]);
  InitFromArrayInLogEkin(np, logEkin.data(), response, opt, valbeg, valend);
}

FCSReturnCode
TFCSEnergyInterpolationSpline::simulate(TFCSSimulationState &simulstate,
                                        const TFCSTruthState *truth,
                                        const TFCSExtrapolationState *) const {
  float Emean;
  float Einit;
  const float Ekin = truth->Ekin();
  if (OnlyScaleEnergy())
    Einit = simulstate.E();
  else
    Einit = Ekin;
  // catch very small values of Ekin (use 1 keV here) and fix the spline lookup
  // to the 1keV value
  const float logEkin =
      (Ekin > Gaudi::Units::keV ? TMath::Log(Ekin)
                                : TMath::Log(Gaudi::Units::keV));
  if (logEkin < m_spline.GetXmin()) {
    Emean = m_spline.Eval(m_spline.GetXmin()) * Einit;
  } else {
    if (logEkin > m_spline.GetXmax()) {
      Emean = (m_spline.Eval(m_spline.GetXmax()) +
               m_spline.Derivative(m_spline.GetXmax()) *
                   (logEkin - m_spline.GetXmax())) *
              Einit;
    } else {
      Emean = m_spline.Eval(logEkin) * Einit;
    }
  }

  if (OnlyScaleEnergy()) {
    ATH_MSG_DEBUG("set E=" << Emean << " for true Ekin=" << truth->Ekin()
                           << " and E=" << Einit);
  } else {
    ATH_MSG_DEBUG("set E=" << Emean << " for true Ekin=" << truth->Ekin());
  }
  simulstate.set_E(Emean);

  return FCSSuccess;
}

void TFCSEnergyInterpolationSpline::Print(Option_t *option) const {
  TString opt(option);
  bool shortprint = opt.Index("short") >= 0;
  bool longprint = msgLvl(MSG::DEBUG) || (msgLvl(MSG::INFO) && !shortprint);
  TString optprint = opt;
  optprint.ReplaceAll("short", "");
  TFCSParametrization::Print(option);

  if (longprint)
    ATH_MSG_INFO(optprint << (OnlyScaleEnergy() ? "  E()*" : "  Ekin()*")
                          << "Spline N=" << m_spline.GetNp() << " "
                          << m_spline.GetXmin()
                          << "<=log(Ekin)<=" << m_spline.GetXmax() << " "
                          << TMath::Exp(m_spline.GetXmin())
                          << "<=Ekin<=" << TMath::Exp(m_spline.GetXmax()));
}

void TFCSEnergyInterpolationSpline::unit_test(
    TFCSSimulationState *simulstate, TFCSTruthState *truth,
    const TFCSExtrapolationState *extrapol, TGraph *grspline) {
  if (!simulstate)
    simulstate = new TFCSSimulationState();
  if (!truth)
    truth = new TFCSTruthState();
  if (!extrapol)
    extrapol = new TFCSExtrapolationState();

  if (!grspline) {
    const int Graph0_n = 9;
    Double_t Graph0_fx1001[Graph0_n] = {1.024,  2.048,  4.094,   8.192,  16.384,
                                        32.768, 65.536, 131.072, 262.144};
    for (int i = 0; i < Graph0_n; ++i)
      Graph0_fx1001[i] *= 1000;

    Double_t Graph0_fy1001[Graph0_n] = {0.6535402, 0.6571529, 0.6843001,
                                        0.7172835, 0.7708416, 0.798819,
                                        0.8187628, 0.8332745, 0.8443931};
    grspline = new TGraph(Graph0_n, Graph0_fx1001, Graph0_fy1001);
  }

  TGraph *grdraw = (TGraph *)grspline->Clone();
  grdraw->SetMarkerColor(46);
  grdraw->SetMarkerStyle(8);

  TFCSEnergyInterpolationSpline test("testTFCSEnergyInterpolationSpline",
                                     "test TFCSEnergyInterpolationSpline");
  test.set_pdgid(22);
  test.set_Ekin_nominal(
      0.5 * (grdraw->GetX()[0] + grdraw->GetX()[grdraw->GetN() - 1]));
  test.set_Ekin_min(grdraw->GetX()[0]);
  test.set_Ekin_max(grdraw->GetX()[grdraw->GetN() - 1]);
  test.set_eta_nominal(0.225);
  test.set_eta_min(0.2);
  test.set_eta_max(0.25);
  test.InitFromArrayInEkin(grspline->GetN(), grspline->GetX(), grspline->GetY(),
                           "b2e2", 0, 0);
  // test.set_OnlyScaleEnergy();
  test.Print();

  truth->set_pdgid(22);

  TGraph *gr = new TGraph();
  gr->SetNameTitle("testTFCSEnergyInterpolationSplineLogX",
                   "test TFCSEnergyInterpolationSpline log x-axis");
  gr->GetXaxis()->SetTitle("Ekin [MeV]");
  gr->GetYaxis()->SetTitle("<E(reco)>/Ekin(true)");

  int ip = 0;
  for (float Ekin = test.Ekin_min() * 0.25; Ekin <= test.Ekin_max() * 4;
       Ekin *= 1.05) {
    // Init LorentzVector for truth. For photon Ekin=E
    truth->SetPxPyPzE(Ekin, 0, 0, Ekin);
    simulstate->set_E(Ekin);
    if (test.simulate(*simulstate, truth, extrapol) != FCSSuccess) {
      return;
    }
    gr->SetPoint(ip, Ekin, simulstate->E() / Ekin);
    ++ip;
  }

  TCanvas *c = new TCanvas(gr->GetName(), gr->GetTitle());
  gr->Draw("APL");
  grdraw->Draw("Psame");
  c->SetLogx();

}
