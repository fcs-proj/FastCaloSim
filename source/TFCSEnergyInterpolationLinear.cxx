/*
  Copyright (C) 2002-2019 CERN for the benefit of the ATLAS collaboration
*/

#include <iostream>

#include "FastCaloSim/TFCSEnergyInterpolationLinear.h"

#include "FastCaloSim/TFCSExtrapolationState.h"
#include "FastCaloSim/TFCSSimulationState.h"
#include "FastCaloSim/TFCSTruthState.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TGraph.h"

//=============================================
//======= TFCSEnergyInterpolationLinear =========
//=============================================

TFCSEnergyInterpolationLinear::TFCSEnergyInterpolationLinear(const char* name,
                                                             const char* title)
    : TFCSParametrization(name, title)
    , m_slope(1)
    , m_offset(0)
{
}

FCSReturnCode TFCSEnergyInterpolationLinear::simulate(
    TFCSSimulationState& simulstate,
    const TFCSTruthState* truth,
    const TFCSExtrapolationState*) const
{
  const float Emean = m_slope * truth->Ekin() + m_offset;

  ATH_MSG_DEBUG("set E=" << Emean << " for true Ekin=" << truth->Ekin());
  simulstate.set_E(Emean);

  return FCSSuccess;
}

void TFCSEnergyInterpolationLinear::Print(Option_t* option) const
{
  TString opt(option);
  bool shortprint = opt.Index("short") >= 0;
  bool longprint = msgLvl(MSG::DEBUG) || (msgLvl(MSG::INFO) && !shortprint);
  TString optprint = opt;
  optprint.ReplaceAll("short", "");
  TFCSParametrization::Print(option);

  if (longprint)
    ATH_MSG_INFO(optprint << "  Emean=" << m_slope << "*Ekin(true) + "
                          << m_offset);
}

void TFCSEnergyInterpolationLinear::unit_test(
    TFCSSimulationState* simulstate,
    TFCSTruthState* truth,
    const TFCSExtrapolationState* extrapol)
{
  if (!simulstate)
    simulstate = new TFCSSimulationState();
  if (!truth)
    truth = new TFCSTruthState();
  if (!extrapol)
    extrapol = new TFCSExtrapolationState();

  TFCSEnergyInterpolationLinear test("testTFCSEnergyInterpolation",
                                     "test linear TFCSEnergyInterpolation");
  test.set_pdgid(22);
  test.set_Ekin_nominal(1000);
  test.set_Ekin_min(1000);
  test.set_Ekin_max(100000);
  test.set_eta_nominal(0.225);
  test.set_eta_min(0.2);
  test.set_eta_max(0.25);
  test.set_slope(0.95);
  test.set_offset(-50);
  test.Print();

  truth->set_pdgid(22);

  TGraph* gr = new TGraph();
  gr->SetNameTitle("testTFCSEnergyInterpolation",
                   "test linear TFCSEnergyInterpolation");
  gr->GetXaxis()->SetTitle("Ekin [MeV]");
  gr->GetYaxis()->SetTitle("<E(reco)>/Ekin(true)");

  int ip = 0;
  for (float Ekin = 1000; Ekin <= 100000; Ekin *= 1.1) {
    // Init LorentzVector for truth. For photon Ekin=E
    truth->SetPxPyPzE(Ekin, 0, 0, Ekin);
    if (test.simulate(*simulstate, truth, extrapol) != FCSSuccess) {
      return;
    }
    gr->SetPoint(ip, Ekin, simulstate->E() / Ekin);
    ++ip;
  }

  TCanvas* c = new TCanvas("testTFCSEnergyInterpolation",
                           "test linear TFCSEnergyInterpolation");
  gr->Draw("APL");
  c->SetLogx();
}
