// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

///////////////////////////////////////////////////////////////////
// TFCSGANEtaSlice.cxx, (c) ATLAS Detector software             //
///////////////////////////////////////////////////////////////////

// class header include
#include <charconv>
#include <cmath>
#include <cstddef>
#include <format>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <system_error>
#include <vector>

#include "FastCaloSim/Core/TFCSGANEtaSlice.h"

#include "CLHEP/Random/RandGauss.h"
#include "FastCaloSim/Core/TFCSNetworkFactory.h"
#include "TFile.h"
#include "TFitResult.h"
#include "TH1D.h"
#include "TTree.h"

TFCSGANEtaSlice::TFCSGANEtaSlice() = default;

// NOLINTBEGIN(modernize-pass-by-value): keep the public signature
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
// NOLINTEND(modernize-pass-by-value)

TFCSGANEtaSlice::~TFCSGANEtaSlice()
{
  // Deleting a nullptr is a noop
  delete m_gan_all;
  delete m_gan_low;
  delete m_gan_high;
}

auto TFCSGANEtaSlice::GetNetAll() const -> VNetworkBase*
{
  if (m_net_all != nullptr) {
    return m_net_all.get();
  }
  return m_gan_all;
}
auto TFCSGANEtaSlice::GetNetLow() const -> VNetworkBase*
{
  if (m_net_low != nullptr) {
    return m_net_low.get();
  }
  return m_gan_low;
}
auto TFCSGANEtaSlice::GetNetHigh() const -> VNetworkBase*
{
  if (m_net_high != nullptr) {
    return m_net_high.get();
  }
  return m_gan_high;
}

auto TFCSGANEtaSlice::IsGanCorrectlyLoaded() const -> bool
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

auto TFCSGANEtaSlice::LoadGAN() -> bool
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
    FCS_MSG_DEBUG("Gan input file name " << inputFileName);
    m_net_all = TFCSNetworkFactory::create(std::move(inputFileName));
    if (m_net_all == nullptr) {
      success = false;
    }
  } else if (m_pid == 2212) {
    inputFileName = m_param.GetInputFolder() + "/neural_net_"
        + std::to_string(m_pid) + "_eta_" + std::to_string(m_etaMin) + "_"
        + std::to_string(m_etaMax) + "_High10.*";
    FCS_MSG_DEBUG("Gan input file name " << inputFileName);
    m_net_all = TFCSNetworkFactory::create(std::move(inputFileName));
    if (m_net_all == nullptr) {
      success = false;
    }
  } else {
    inputFileName = m_param.GetInputFolder() + "/neural_net_"
        + std::to_string(m_pid) + "_eta_" + std::to_string(m_etaMin) + "_"
        + std::to_string(m_etaMax) + "_High12.*";
    FCS_MSG_DEBUG("Gan input file name " << inputFileName);
    m_net_high = TFCSNetworkFactory::create(std::move(inputFileName));
    if (m_net_high == nullptr) {
      success = false;
    }

    inputFileName = m_param.GetInputFolder() + "/neural_net_"
        + std::to_string(m_pid) + "_eta_" + std::to_string(m_etaMin) + "_"
        + std::to_string(m_etaMax) + "_UltraLow12.*";
    m_net_low = TFCSNetworkFactory::create(std::move(inputFileName));
    if (m_net_low == nullptr) {
      success = false;
    }
  }
  return success;
}

void TFCSGANEtaSlice::CalculateMeanPointFromDistributionOfR()
{
  std::string const rootFileName = m_param.GetInputFolder() + "/rootFiles/pid"
      + std::to_string(m_pid) + "_E1048576_eta_" + std::to_string(m_etaMin)
      + "_" + std::to_string(m_etaMin + 5) + ".root";
  FCS_MSG_DEBUG("Opening file " << rootFileName);
  std::unique_ptr<TFile> file(TFile::Open(rootFileName.c_str(), "read"));
  if (!file || file->IsZombie()) {
    throw std::runtime_error("Failed to open or initialize ROOT file: "
                             + rootFileName);
  }
  for (int const layer : m_param.GetRelevantLayers()) {
    FCS_MSG_DEBUG("Layer " << layer);
    TFCSGANXMLParameters::Binning binsInLayers = m_param.GetBinning();
    TH2D* h2 = &binsInLayers[layer];

    std::string histoName = "r" + std::to_string(layer) + "w";
    TH1D* h1 = static_cast<TH1D*>(file->Get(histoName.c_str()));
    if (h1 == nullptr || std::isnan(h1->Integral())) {
      histoName = "r" + std::to_string(layer);
      h1 = static_cast<TH1D*>(file->Get(histoName.c_str()));
    }
    if (h1 == nullptr) {
      throw std::runtime_error(
          std::format("Failed to load histogram {} from ROOT file: {}",
                      histoName,
                      rootFileName));
    }

    TAxis const* x = h2->GetXaxis();
    for (int ix = 1; ix <= h2->GetNbinsX(); ++ix) {
      FCS_MSG_DEBUG(ix);
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
  FCS_MSG_DEBUG("Done initializing fits");
}

void TFCSGANEtaSlice::ExtractExtrapolatorMeansFromInputs()
{
  std::string const rootFileName = m_param.GetInputFolder() + "/rootFiles/pid"
      + std::to_string(m_pid) + "_E65536_eta_" + std::to_string(m_etaMin) + "_"
      + std::to_string(m_etaMin + 5) + "_validation.root";
  FCS_MSG_DEBUG("Opening file " << rootFileName);
  std::unique_ptr<TFile> file(TFile::Open(rootFileName.c_str(), "read"));
  if (!file || file->IsZombie()) {
    throw std::runtime_error("Failed to open or initialize ROOT file: "
                             + rootFileName);
  }
  auto* tree = file->Get<TTree>("rootTree");
  if (tree == nullptr) {
    throw std::runtime_error("No rootTree in " + rootFileName);
  }
  for (int const layer : m_param.GetRelevantLayers()) {
    std::string const branchName = "extrapWeight_" + std::to_string(layer);
    // TTree::Draw fills the histogram that it finds by name in the current
    // directory (the file). Deleting it removes it from there again.
    auto h = std::make_unique<TH1D>("h", "h", 100, 0.01, 1);
    std::string const command = branchName + ">>h";
    tree->Draw(command.c_str());
    m_extrapolatorWeights[layer] = h->GetMean();
    FCS_MSG_DEBUG("Extrapolation: layer " << layer << " mean "
                                          << m_extrapolatorWeights[layer]);
  }
}

auto TFCSGANEtaSlice::GetVoxelEnergyFractions(
    const TFCSTruthState* truth,
    const TFCSExtrapolationState* /*extrapol*/,
    TFCSSimulationState& simulstate) const -> std::vector<double>
{
  double randUniformZ = 0.;
  NetworkInputs inputs;

  int maxExp = 0;
  int minExp = 0;
  if (m_pid == 22 || std::abs(m_pid) == 11) {
    if (truth->P() > 4096)
    {  // This is the momentum, not the energy, because the split is
       // based on the samples which are produced with the momentum
      maxExp = 22;
      minExp = 12;
    } else {
      maxExp = 12;
      minExp = 6;
    }
  } else if (std::abs(m_pid) == 211) {
    maxExp = 22;
    minExp = 8;
  } else if (std::abs(m_pid) == 2212) {
    maxExp = 22;
    minExp = 10;
  }

  const auto p_min = static_cast<int>(std::pow(2, minExp));
  const auto p_max = static_cast<int>(std::pow(2, maxExp));
  // Keep min and max without mass offset as we do not train on antiparticles
  double const Ekin_min =
      std::sqrt(std::pow(p_min, 2) + std::pow(truth->M(), 2)) - truth->M();
  double const Ekin_max =
      std::sqrt(std::pow(p_max, 2) + std::pow(truth->M(), 2)) - truth->M();

  for (int i = 0; i < m_param.GetLatentSpaceSize(); i++) {
    randUniformZ = CLHEP::RandGauss::shoot(simulstate.randomEngine(), 0.5, 0.5);
    inputs["Noise"].insert(std::pair<std::string, double>(
        "variable_" + std::to_string(i), randUniformZ));
  }

  inputs["mycond"].insert(std::pair<std::string, double>(
      "variable_0", log(truth->Ekin() / Ekin_min) / log(Ekin_max / Ekin_min)));

  if (m_param.GetGANVersion() >= 2) {
    // Conditioning on |eta| would only be needed in transition regions. All
    // current GANs have this input, so it is filled with zero.
    inputs["mycond"].insert(std::pair<std::string, double>("variable_1", 0));
  }

  VNetworkBase::NetworkOutputs outputs;
  if (m_param.GetGANVersion() == 1 || m_pid == 211 || m_pid == 2212) {
    outputs = GetNetAll()->compute(inputs);
  } else {
    if (truth->P() > 4096)
    {  // This is the momentum, not the energy, because the split is
       // based on the samples which are produced with the momentum
      FCS_MSG_DEBUG("Computing outputs given inputs for high");
      outputs = GetNetHigh()->compute(inputs);
    } else {
      outputs = GetNetLow()->compute(inputs);
    }
  }
  FCS_MSG_DEBUG("Start Network inputs ~~~~~~~~");
  FCS_MSG_DEBUG(VNetworkBase::representNetworkInputs(inputs, 10000));
  FCS_MSG_DEBUG("End Network inputs ~~~~~~~~");
  FCS_MSG_DEBUG("Start Network outputs ~~~~~~~~");
  FCS_MSG_DEBUG(VNetworkBase::representNetworkOutputs(outputs, 10000));
  FCS_MSG_DEBUG("End Network outputs ~~~~~~~~");

  // The network names its outputs by the voxel number: index them once
  std::vector<double> fractions(outputs.size());
  for (const auto& [name, value] : outputs) {
    std::size_t voxel = 0;
    const char* const end = name.data() + name.size();
    const auto [ptr, ec] = std::from_chars(name.data(), end, voxel);
    if (ec != std::errc {} || ptr != end || voxel >= fractions.size()) {
      FCS_MSG_ERROR("Network output " << name << " is not a voxel number below "
                                      << fractions.size());
      return {};
    }
    fractions[voxel] = value;
  }
  return fractions;
}

void TFCSGANEtaSlice::Print() const
{
  FCS_MSG_INFO("LWTNN Handler parameters");
  FCS_MSG_INFO("  pid: " << m_pid);
  FCS_MSG_INFO("  etaMin:" << m_etaMin);
  FCS_MSG_INFO("  etaMax: " << m_etaMax);
  m_param.Print();
}

void TFCSGANEtaSlice::fixHists()
{
  m_param.fixHists();
}
