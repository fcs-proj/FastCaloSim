/*
  Copyright (C) 2002-2020 CERN for the benefit of the ATLAS collaboration
*/

#include <iostream>
#include <vector>

#include "FastCaloSim/Core/TFCSEnergyInterpolationHistogram.h"

#include "FastCaloSim/Core/TFCSExtrapolationState.h"
#include "FastCaloSim/Core/TFCSSimulationState.h"
#include "FastCaloSim/Core/TFCSTruthState.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"

namespace Gaudi
{
namespace Units
{
constexpr double megaelectronvolt = 1.;
constexpr double kiloelectronvolt = 1.e-3 * megaelectronvolt;
constexpr double keV = kiloelectronvolt;
}  // namespace Units
}  // namespace Gaudi

//=============================================
//======= TFCSEnergyInterpolation =========
//=============================================

TFCSEnergyInterpolationHistogram::TFCSEnergyInterpolationHistogram(
    const char* name, const char* title)
    : TFCSParametrization(name, title)
{
}

FCSReturnCode TFCSEnergyInterpolationHistogram::simulate(
    TFCSSimulationState& simulstate,
    const TFCSTruthState* truth,
    const TFCSExtrapolationState*) const
{
  float Emean;
  float Einit;
  const float Ekin = truth->Ekin();
  if (OnlyScaleEnergy())
    Einit = simulstate.E();
  else
    Einit = Ekin;
  if (Einit < m_hist.GetXaxis()->GetBinLowEdge(1)) {
    Emean = m_hist.GetBinContent(1) * Einit;
  } else {
    if (Einit > m_hist.GetXaxis()->GetBinUpEdge(m_hist.GetNbinsX())) {
      Emean = m_hist.GetBinContent(m_hist.GetNbinsX()) * Einit;
    } else {
      Emean = m_hist.GetBinContent(m_hist.GetXaxis()->FindBin(Einit)) * Einit;
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

void TFCSEnergyInterpolationHistogram::Print(Option_t* option) const
{
  TString opt(option);
  bool shortprint = opt.Index("short") >= 0;
  bool longprint = msgLvl(MSG::DEBUG) || (msgLvl(MSG::INFO) && !shortprint);
  TString optprint = opt;
  optprint.ReplaceAll("short", "");
  TFCSParametrization::Print(option);

  if (longprint)
    ATH_MSG_INFO(
        optprint << (OnlyScaleEnergy() ? "  E()*" : "  Ekin()*")
                 << "histNbins=" << m_hist.GetNbinsX() << " "
                 << m_hist.GetXaxis()->GetBinLowEdge(1) << "<=Ekin<="
                 << m_hist.GetXaxis()->GetBinUpEdge(m_hist.GetNbinsX()));
}

void TFCSEnergyInterpolationHistogram::unit_test(
    TFCSSimulationState* simulstate,
    TFCSTruthState* truth,
    const TFCSExtrapolationState* extrapol,
    TH1F* hist)
{
  if (!simulstate)
    simulstate = new TFCSSimulationState();
  if (!truth)
    truth = new TFCSTruthState();
  if (!extrapol)
    extrapol = new TFCSExtrapolationState();

  if (!hist) {
    hist = new TH1F("h1", "h1 title", 16, 0., 512.);
    hist->SetBinContent(1, 0.687165);
    hist->SetBinContent(2, 0.756837);
    hist->SetBinContent(3, 0.836673);
    hist->SetBinContent(4, 0.896336);
    hist->SetBinContent(5, 0.889785);
    hist->SetBinContent(6, 0.901266);
    hist->SetBinContent(7, 0.888937);
    hist->SetBinContent(8, 0.919943);
    hist->SetBinContent(9, 0.941806);
    hist->SetBinContent(10, 0.934668);
    hist->SetBinContent(11, 0.939502);
    hist->SetBinContent(12, 0.940796);
    hist->SetBinContent(13, 0.945894);
    hist->SetBinContent(14, 0.955445);
    hist->SetBinContent(15, 0.955593);
    hist->SetBinContent(16, 0.943673);
  }

  TH1F* hdraw = (TH1F*)hist->Clone();
  hdraw->SetMarkerColor(46);
  hdraw->SetMarkerStyle(8);

  TFCSEnergyInterpolationHistogram test(
      "testTFCSEnergyInterpolationHistogram",
      "test TFCSEnergyInterpolationHistogram");
  test.set_pdgid(22);
  test.set_Ekin_nominal(
      0.5
      * (hdraw->GetXaxis()->GetBinLowEdge(1)
         + hdraw->GetXaxis()->GetBinUpEdge(hdraw->GetNbinsX())));
  test.set_Ekin_min(hdraw->GetXaxis()->GetBinLowEdge(1));
  test.set_Ekin_max(hdraw->GetXaxis()->GetBinUpEdge(hdraw->GetNbinsX()));
  test.set_eta_nominal(0.225);
  test.set_eta_min(0.2);
  test.set_eta_max(0.25);
  test.InitFromHist(*hist);
  // test.set_OnlyScaleEnergy();
  test.Print();
  test.hist().Dump();

  truth->set_pdgid(22);

  TGraph* gr = new TGraph();
  gr->SetNameTitle("testTFCSEnergyInterpolationHistogramLogX",
                   "test TFCSEnergyInterpolationHistogram");
  gr->GetXaxis()->SetTitle("Ekin [MeV]");
  gr->GetYaxis()->SetTitle("<E(reco)>/Ekin(true)");

  int ip = 0;
  for (float Ekin = std::max(test.Ekin_min() * 0.25, 0.1);
       Ekin <= test.Ekin_max() * 4;
       Ekin *= 1.05)
  {
    // Init LorentzVector for truth. For photon Ekin=E
    truth->SetPxPyPzE(Ekin, 0, 0, Ekin);
    simulstate->set_E(Ekin);
    if (test.simulate(*simulstate, truth, extrapol) != FCSSuccess) {
      return;
    }
    gr->SetPoint(ip, Ekin, simulstate->E() / Ekin);
    ++ip;
  }

  TCanvas* c = new TCanvas(hdraw->GetName(), hdraw->GetTitle());
  hdraw->Draw("HIST");
  gr->Draw("same,APL");
  c->SetLogx();
}
