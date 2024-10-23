// Copyright (c) 2024 CERN for the benefit of the FastCaloSim project

///////////////////////////////////////////////////////////////////
// TFCSGANEtaSlice.cxx, (c) ATLAS Detector software             //
///////////////////////////////////////////////////////////////////

// class header include
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include "FastCaloSim/Core/TFCSGANEtaSlice.h"

#include "CLHEP/Random/RandGauss.h"
#include "FastCaloSim/Core/TFCSNetworkFactory.h"
#include "TFile.h"
#include "TFitResult.h"
#include "TH1D.h"
#include "TMath.h"
#include "TTree.h"

TFCSGANEtaSlice::TFCSGANEtaSlice() {}

TFCSGANEtaSlice::TFCSGANEtaSlice(int pid,
                                 int etaMin,
                                 int etaMax,
                                 const TFCSGANXMLParameters& param)
    : m_pid(pid)
    , m_etaMin(etaMin)
    , m_etaMax(etaMax)
    , m_param(param)
{
}

TFCSGANEtaSlice::~TFCSGANEtaSlice()
{
  // Deleting a nullptr is a noop
  delete m_gan_all;
  delete m_gan_low;
  delete m_gan_high;
}

VNetworkBase* TFCSGANEtaSlice::GetNetAll() const
{
  if (m_net_all != nullptr)
    return m_net_all.get();
  return m_gan_all;
}
VNetworkBase* TFCSGANEtaSlice::GetNetLow() const
{
  if (m_net_low != nullptr)
    return m_net_low.get();
  return m_gan_low;
}
VNetworkBase* TFCSGANEtaSlice::GetNetHigh() const
{
  if (m_net_high != nullptr)
    return m_net_high.get();
  return m_gan_high;
}

bool TFCSGANEtaSlice::IsGanCorrectlyLoaded() const
{
  if (m_pid == 211 || m_pid == 2212) {
    if (GetNetAll() == nullptr) {
      return false;
    }
  } else {
    if (GetNetHigh() == nullptr || GetNetLow() == nullptr) {
      return false;
    }
  }
  return true;
}

bool TFCSGANEtaSlice::LoadGAN()
{
  // Now load new data
  std::string inputFileName;

  CalculateMeanPointFromDistributionOfR();
  ExtractExtrapolatorMeansFromInputs();

  bool success = true;

  if (m_pid == 211) {
    inputFileName = m_param.GetInputFolder() + "/neural_net_"
        + std::to_string(m_pid) + "_eta_" + std::to_string(m_etaMin) + "_"
        + std::to_string(m_etaMax) + "_All.*";
    ATH_MSG_DEBUG("Gan input file name " << inputFileName);
    m_net_all = TFCSNetworkFactory::create(inputFileName);
    if (m_net_all == nullptr)
      success = false;
  } else if (m_pid == 2212) {
    inputFileName = m_param.GetInputFolder() + "/neural_net_"
        + std::to_string(m_pid) + "_eta_" + std::to_string(m_etaMin) + "_"
        + std::to_string(m_etaMax) + "_High10.*";
    ATH_MSG_DEBUG("Gan input file name " << inputFileName);
    m_net_all = TFCSNetworkFactory::create(inputFileName);
    if (m_net_all == nullptr)
      success = false;
  } else {
    inputFileName = m_param.GetInputFolder() + "/neural_net_"
        + std::to_string(m_pid) + "_eta_" + std::to_string(m_etaMin) + "_"
        + std::to_string(m_etaMax) + "_High12.*";
    ATH_MSG_DEBUG("Gan input file name " << inputFileName);
    m_net_high = TFCSNetworkFactory::create(inputFileName);
    if (m_net_high == nullptr)
      success = false;

    inputFileName = m_param.GetInputFolder() + "/neural_net_"
        + std::to_string(m_pid) + "_eta_" + std::to_string(m_etaMin) + "_"
        + std::to_string(m_etaMax) + "_UltraLow12.*";
    m_net_low = TFCSNetworkFactory::create(inputFileName);
    if (m_net_low == nullptr)
      success = false;
  }
  return success;
}

void TFCSGANEtaSlice::CalculateMeanPointFromDistributionOfR()
{
  std::string rootFileName = m_param.GetInputFolder() + "/rootFiles/pid"
      + std::to_string(m_pid) + "_E1048576_eta_" + std::to_string(m_etaMin)
      + "_" + std::to_string(m_etaMin + 5) + ".root";
  ATH_MSG_DEBUG("Opening file " << rootFileName);
  TFile* file = TFile::Open(rootFileName.c_str(), "read");
  for (int layer : m_param.GetRelevantLayers()) {
    ATH_MSG_DEBUG("Layer " << layer);
    TFCSGANXMLParameters::Binning binsInLayers = m_param.GetBinning();
    TH2D* h2 = &binsInLayers[layer];

    std::string histoName = "r" + std::to_string(layer) + "w";
    TH1D* h1 = (TH1D*)file->Get(histoName.c_str());
    if (TMath::IsNaN(h1->Integral())) {
      histoName = "r" + std::to_string(layer);
      h1 = (TH1D*)file->Get(histoName.c_str());
    }

    TAxis* x = (TAxis*)h2->GetXaxis();
    for (int ix = 1; ix <= h2->GetNbinsX(); ++ix) {
      ATH_MSG_DEBUG(ix);
      h1->GetXaxis()->SetRangeUser(x->GetBinLowEdge(ix), x->GetBinUpEdge(ix));

      double result = 0;
      if (h1->Integral() > 0 && h1->GetNbinsX() > 2) {
        TFitResultPtr res(0);

        res = h1->Fit("expo", "SQ");
        if (res >= 0 && !std::isnan(res->Parameter(0))) {
          result = res->Parameter(1);
        }
      }
      m_allFitResults[layer].push_back(result);
    }
  }
  ATH_MSG_DEBUG("Done initializing fits");
}

void TFCSGANEtaSlice::ExtractExtrapolatorMeansFromInputs()
{
  std::string rootFileName = m_param.GetInputFolder() + "/rootFiles/pid"
      + std::to_string(m_pid) + "_E65536_eta_" + std::to_string(m_etaMin) + "_"
      + std::to_string(m_etaMin + 5) + "_validation.root";
  ATH_MSG_DEBUG("Opening file " << rootFileName);
  TFile* file = TFile::Open(rootFileName.c_str(), "read");
  for (int layer : m_param.GetRelevantLayers()) {
    std::string branchName = "extrapWeight_" + std::to_string(layer);
    TH1D* h = new TH1D("h", "h", 100, 0.01, 1);
    TTree* tree = (TTree*)file->Get("rootTree");
    std::string command = branchName + ">>h";
    tree->Draw(command.c_str());
    m_extrapolatorWeights[layer] = h->GetMean();
    ATH_MSG_DEBUG("Extrapolation: layer " << layer << " mean "
                                          << m_extrapolatorWeights[layer]);
  }
}

VNetworkBase::NetworkOutputs TFCSGANEtaSlice::GetNetworkOutputs(
    const TFCSTruthState* truth,
    const TFCSExtrapolationState* extrapol,
    TFCSSimulationState simulstate) const
{
  double randUniformZ = 0.;
  NetworkInputs inputs;

  int maxExp = 0, minExp = 0;
  if (m_pid == 22 || fabs(m_pid) == 11) {
    if (truth->P() > 4096)
    {  // This is the momentum, not the energy, because the split is
       // based on the samples which are produced with the momentum
      maxExp = 22;
      minExp = 12;
    } else {
      maxExp = 12;
      minExp = 6;
    }
  } else if (fabs(m_pid) == 211) {
    maxExp = 22;
    minExp = 8;
  } else if (fabs(m_pid) == 2212) {
    maxExp = 22;
    minExp = 10;
  }

  int p_min = std::pow(2, minExp);
  int p_max = std::pow(2, maxExp);
  // Keep min and max without mass offset as we do not train on antiparticles
  double Ekin_min =
      std::sqrt(std::pow(p_min, 2) + std::pow(truth->M(), 2)) - truth->M();
  double Ekin_max =
      std::sqrt(std::pow(p_max, 2) + std::pow(truth->M(), 2)) - truth->M();

  for (int i = 0; i < m_param.GetLatentSpaceSize(); i++) {
    randUniformZ = CLHEP::RandGauss::shoot(simulstate.randomEngine(), 0.5, 0.5);
    inputs["Noise"].insert(std::pair<std::string, double>(
        "variable_" + std::to_string(i), randUniformZ));
  }

  // double e = log(truth->Ekin()/Ekin_min)/log(Ekin_max/Ekin_min) ;
  // Could be uncommented , but would need the line above too
  // ATH_MSG_DEBUG( "Check label: " << e <<" Ekin:" << truth->Ekin() <<" p:" <<
  //                truth->P() <<" mass:" << truth->M() <<" Ekin_off:" <<
  //                truth->Ekin_off() << " Ekin_min:"<<Ekin_min<<"
  //                Ekin_max:"<<Ekin_max);
  // inputs["mycond"].insert ( std::pair<std::string,double>("variable_0",
  // truth->Ekin()/(std::pow(2,maxExp))) ); //Old conditioning using linear
  // interpolation, now use logaritminc interpolation
  inputs["mycond"].insert(std::pair<std::string, double>(
      "variable_0", log(truth->Ekin() / Ekin_min) / log(Ekin_max / Ekin_min)));

  if (m_param.GetGANVersion() >= 2) {
    if (false) {  // conditioning on eta, should only be needed in transition
                  // regions and added only to the GANs that use it, for now all
                  // GANs have 3 conditioning inputs so filling zeros
      inputs["mycond"].insert(std::pair<std::string, double>(
          "variable_1", fabs(extrapol->IDCaloBoundary_eta())));
    } else {
      inputs["mycond"].insert(std::pair<std::string, double>("variable_1", 0));
    }
  }

  VNetworkBase::NetworkOutputs outputs;
  if (m_param.GetGANVersion() == 1 || m_pid == 211 || m_pid == 2212) {
    outputs = GetNetAll()->compute(inputs);
  } else {
    if (truth->P() > 4096)
    {  // This is the momentum, not the energy, because the split is
       // based on the samples which are produced with the momentum
      ATH_MSG_DEBUG("Computing outputs given inputs for high");
      outputs = GetNetHigh()->compute(inputs);
    } else {
      outputs = GetNetLow()->compute(inputs);
    }
  }
  ATH_MSG_DEBUG("Start Network inputs ~~~~~~~~");
  ATH_MSG_DEBUG(VNetworkBase::representNetworkInputs(inputs, 10000));
  ATH_MSG_DEBUG("End Network inputs ~~~~~~~~");
  ATH_MSG_DEBUG("Start Network outputs ~~~~~~~~");
  ATH_MSG_DEBUG(VNetworkBase::representNetworkOutputs(outputs, 10000));
  ATH_MSG_DEBUG("End Network outputs ~~~~~~~~");
  return outputs;
}

void TFCSGANEtaSlice::Print() const
{
  ATH_MSG_INFO("LWTNN Handler parameters");
  ATH_MSG_INFO("  pid: " << m_pid);
  ATH_MSG_INFO("  etaMin:" << m_etaMin);
  ATH_MSG_INFO("  etaMax: " << m_etaMax);
  m_param.Print();
}
