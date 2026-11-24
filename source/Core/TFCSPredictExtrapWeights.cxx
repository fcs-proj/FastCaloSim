// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <iostream>

#include "FastCaloSim/Core/TFCSPredictExtrapWeights.h"

#include <CLHEP/Random/RanluxEngine.h>

#include "FastCaloSim/Core/TFCSExtrapolationState.h"
#include "FastCaloSim/Core/TFCSSimulationState.h"
#include "FastCaloSim/Core/TFCSTruthState.h"
#include "FastCaloSim/Definitions/Index.h"
#include "FastCaloSim/Definitions/Phi.h"
#include "FastCaloSim/Geometry/CaloGeo.h"
#include "TBuffer.h"
#include "TClass.h"
// LWTNN
#include "lwtnn/LightweightNeuralNetwork.hh"
#include "lwtnn/parse_json.hh"

// XML reader
#include <libxml/parser.h>
#include <libxml/tree.h>
#include <libxml/xmlmemory.h>
#include <libxml/xmlreader.h>
#include <libxml/xpath.h>
#include <libxml/xpathInternals.h>

//=============================================
//======= TFCSPredictExtrapWeights =========
//=============================================

namespace
{
// AuxInfo index of the predicted extrapolation weight of a layer
auto extrap_weight_index(std::size_t layer) -> std::uint32_t
{
  return TFCSSimulationState::layerAuxIndex("ExtrapWeight"_FCShash, layer);
}
}  // namespace

TFCSPredictExtrapWeights::TFCSPredictExtrapWeights(const char* name,
                                                   const char* title,
                                                   CaloGeo* geo)
    : TFCSLateralShapeParametrizationHitBase(name, title)
    , m_geo(geo)
{
  set_freemem();
  set_UseHardcodedWeight();
}

// Destructor
TFCSPredictExtrapWeights::~TFCSPredictExtrapWeights()
{
  delete m_input;

  delete m_relevantLayers;

  delete m_normLayers;

  delete m_normMeans;

  delete m_normStdDevs;
}

auto TFCSPredictExtrapWeights::operator==(
    const TFCSParametrizationBase& ref) const -> bool
{
  if (IsA() != ref.IsA()) {
    FCS_MSG_DEBUG("operator==: different class types "
                  << IsA()->GetName() << " != " << ref.IsA()->GetName());
    return false;
  }
  const auto& ref_typed = static_cast<const TFCSPredictExtrapWeights&>(ref);

  if (TFCSParametrizationBase::compare(ref)) {
    return true;
  }
  if (!TFCSParametrization::compare(ref)) {
    return false;
  }
  if (!TFCSLateralShapeParametrization::compare(ref)) {
    return false;
  }

  // Without the network configuration (never loaded, or freed after reading
  // with freemem()) it is unknown whether the networks are the same
  if (m_input == nullptr || ref_typed.m_input == nullptr) {
    return false;
  }
  return *m_input == *ref_typed.m_input;
}

// getNormInputs()
// Get values needed to normalize inputs
auto TFCSPredictExtrapWeights::getNormInputs(
    const std::string& etaBin, const std::string& FastCaloTXTInputFolderName)
    -> bool
{
  FCS_MSG_DEBUG(" Getting normalization inputs... ");

  // Open corresponding TXT file and extract mean/std dev for each variable
  if (m_normLayers != nullptr) {
    m_normLayers->clear();
  } else {
    m_normLayers = new std::vector<int>();
  }
  if (m_normMeans != nullptr) {
    m_normMeans->clear();
  } else {
    m_normMeans = new std::vector<float>();
  }
  if (m_normStdDevs != nullptr) {
    m_normStdDevs->clear();
  } else {
    m_normStdDevs = new std::vector<float>();
  }
  std::string const inputFileName = FastCaloTXTInputFolderName
      + "MeanStdDevEnergyFractions_eta_" + etaBin + ".txt";
  FCS_MSG_DEBUG(" Opening " << inputFileName);
  std::ifstream inputTXT(inputFileName);
  if (inputTXT.is_open()) {
    std::string line;
    while (getline(inputTXT, line)) {
      std::stringstream ss(line);
      unsigned int counter = 0;
      while (ss.good()) {
        std::string substr;
        getline(ss, substr, ' ');
        if (counter == 0) {  // Get index (#layer or -1 if var == etrue)
          if (substr != "etrue") {
            int const index = std::stoi(substr.substr(substr.find('_') + 1));
            m_normLayers->push_back(index);
          } else {  // etrue
            m_normLayers->push_back(-1);
          }
        } else if (counter == 1) {
          m_normMeans->push_back(std::stof(substr));
        } else if (counter == 2) {
          m_normStdDevs->push_back(std::stof(substr));
        }
        counter++;
      }
    }
    inputTXT.close();
  } else {
    FCS_MSG_ERROR(" Unable to open file " << inputFileName);
    return false;
  }

  return true;
}

// prepareInputs()
// Prepare input variables to the Neural Network
auto TFCSPredictExtrapWeights::prepareInputs(TFCSSimulationState& simulstate,
                                             const double truthE) const
    -> std::map<std::string, double>
{
  // The normalization constants are stored in single precision
  const auto normalize = [this](double value, std::size_t index) -> double
  {
    return (value - static_cast<double>(std::as_const(m_normMeans)->at(index)))
        / static_cast<double>(std::as_const(m_normStdDevs)->at(index));
  };

  std::map<std::string, double> inputVariables;
  for (const std::size_t ilayer : m_geo->layer_ids()) {
    if (std::ranges::contains(*m_relevantLayers, ilayer)) {
      const std::string layer = std::to_string(ilayer);
      // Find index
      auto itr = std::ranges::find(std::as_const(*m_normLayers), ilayer);
      if (itr != m_normLayers->cend()) {
        const auto index = FastCaloSim::index_of(*m_normLayers, itr);
        inputVariables["ef_" + layer] =
            normalize(simulstate.Efrac(ilayer), index);
      } else {
        FCS_MSG_ERROR("Normalization information not found for layer "
                      << ilayer);
      }
    }
  }
  // Find index for truth energy
  auto itr = std::ranges::find(std::as_const(*m_normLayers), -1);
  const auto index = FastCaloSim::index_of(*m_normLayers, itr);
  inputVariables["etrue"] = normalize(truthE, index);
  if (is_match_pdgid(22)) {
    inputVariables["pdgId"] = 1;  // one hot encoding
  } else if (is_match_pdgid(11) || is_match_pdgid(-11)) {
    inputVariables["pdgId"] = 0;  // one hot encoding
  }
  return inputVariables;
}

// simulate()
// get predicted extrapolation weights and save them as AuxInfo in simulstate
auto TFCSPredictExtrapWeights::simulate(
    TFCSSimulationState& simulstate,
    const TFCSTruthState* truth,
    const TFCSExtrapolationState* /*extrapol*/) const -> FCSReturnCode
{
  // Get inputs to Neural Network
  std::map<std::string, double> const inputVariables =
      prepareInputs(simulstate, truth->E() * 0.001);

  // Get predicted extrapolation weights
  auto outputs = m_nn->compute(inputVariables);
  for (const std::size_t ilayer : m_geo->layer_ids()) {
    if (std::ranges::contains(*m_relevantLayers, ilayer)) {
      FCS_MSG_DEBUG("TFCSPredictExtrapWeights::simulate: layer: "
                    << ilayer << " weight: "
                    << outputs["extrapWeight_" + std::to_string(ilayer)]);
      double weight = outputs["extrapWeight_" + std::to_string(ilayer)];
      // Protections
      if (weight < 0) {
        weight = 0;
      } else if (weight > 1) {
        weight = 1;
      }
      simulstate.setAuxInfo<double>(extrap_weight_index(ilayer), weight);
    } else {  // use weight=0.5 for non-relevant layers
      FCS_MSG_DEBUG(
          "Setting weight=0.5 for layer = " << std::to_string(ilayer));
      simulstate.setAuxInfo<double>(extrap_weight_index(ilayer), 0.5);
    }
  }
  return FCSSuccess;
}

// simulate_hit()
auto TFCSPredictExtrapWeights::simulate_hit(
    Hit& hit,
    TFCSSimulationState& simulstate,
    const TFCSTruthState* /*truth*/,
    const TFCSExtrapolationState* extrapol) -> FCSReturnCode
{
  const std::size_t cs = calosample();

  // Get corresponding predicted extrapolation weight from simulstate
  if (!simulstate.hasAuxInfo(extrap_weight_index(cs))) {
    // missing AuxInfo
    FCS_MSG_FATAL(
        "Simulstate is not decorated with extrapolation weights for cs = "
        << std::to_string(cs));
    return FCSFatal;
  }
  const auto extrapWeight =
      simulstate.getAuxInfo<double>(extrap_weight_index(cs));

  double eta = (1. - extrapWeight) * extrapol->eta(cs, FastCaloSim::SubPos::ENT)
      + extrapWeight * extrapol->eta(cs, FastCaloSim::SubPos::EXT);
  double phi =
      FastCaloSim::interpolate_phi(extrapol->phi(cs, FastCaloSim::SubPos::ENT),
                                   extrapol->phi(cs, FastCaloSim::SubPos::EXT),
                                   extrapWeight);
  double extrapWeight_for_r_z = extrapWeight;
  if (UseHardcodedWeight()) {
    extrapWeight_for_r_z = 0.5;
    FCS_MSG_DEBUG(
        "Will use extrapWeight=0.5 for r and z when constructing a hit");
  } else {
    FCS_MSG_DEBUG(
        "Will use predicted extrapWeight also for r and z when "
        "constructing a hit");
  }
  double r =
      (1. - extrapWeight_for_r_z) * extrapol->r(cs, FastCaloSim::SubPos::ENT)
      + extrapWeight_for_r_z * extrapol->r(cs, FastCaloSim::SubPos::EXT);
  double z =
      (1. - extrapWeight_for_r_z) * extrapol->z(cs, FastCaloSim::SubPos::ENT)
      + extrapWeight_for_r_z * extrapol->z(cs, FastCaloSim::SubPos::EXT);

  if (!std::isfinite(r) || !std::isfinite(z) || !std::isfinite(eta)
      || !std::isfinite(phi))
  {
    FCS_MSG_WARNING(
        "Extrapolator contains NaN or infinite number.\nSetting "
        "center position to calo boundary.");
    FCS_MSG_WARNING("Before fix: center_r: "
                    << r << " center_z: " << z << " center_phi: " << phi
                    << " center_eta: " << eta << " weight: " << extrapWeight
                    << " cs: " << cs);
    // If extrapolator fails we can set position to calo boundary
    r = extrapol->IDCaloBoundary_r();
    z = extrapol->IDCaloBoundary_z();
    eta = extrapol->IDCaloBoundary_eta();
    phi = extrapol->IDCaloBoundary_phi();
    FCS_MSG_WARNING("After fix: center_r: "
                    << r << " center_z: " << z << " center_phi: " << phi
                    << " center_eta: " << eta << " weight: " << extrapWeight
                    << " cs: " << cs);
  }

  hit.setCenter_r(r);
  hit.setCenter_z(z);
  hit.setCenter_eta(eta);
  hit.setCenter_phi(phi);

  FCS_MSG_DEBUG("TFCSPredictExtrapWeights: center_r: "
                << hit.center_r() << " center_z: " << hit.center_z()
                << " center_phi: " << hit.center_phi()
                << " center_eta: " << hit.center_eta()
                << " weight: " << extrapWeight << " cs: " << cs);

  return FCSSuccess;
}

// initializeNetwork()
// Initialize lwtnn network
auto TFCSPredictExtrapWeights::initializeNetwork(
    int pid,
    const std::string& etaBin,
    const std::string& FastCaloNNInputFolderName) -> bool
{
  FCS_MSG_INFO(
      "Using FastCaloNNInputFolderName: " << FastCaloNNInputFolderName);
  set_pdgid(pid);

  std::string const inputFileName =
      FastCaloNNInputFolderName + "NN_" + etaBin + ".json";
  FCS_MSG_DEBUG("Will read JSON file: " << inputFileName);
  if (inputFileName.empty()) {
    FCS_MSG_ERROR("Could not find json file " << inputFileName);
    return false;
  }
  FCS_MSG_INFO("For pid: " << pid << " and etaBin" << etaBin
                           << ", loading json file " << inputFileName);
  std::ifstream input(inputFileName);
  std::stringstream sin;
  sin << input.rdbuf();
  input.close();
  auto config = lwt::parse_json(sin);
  m_nn = std::make_unique<lwt::LightweightNeuralNetwork>(
      config.inputs, config.layers, config.outputs);

  delete m_input;

  m_input = new std::string(sin.str());
  // Extract relevant layers from the outputs
  m_relevantLayers = new std::vector<int>();
  for (auto name : config.outputs) {
    int const layer = std::stoi(
        name.erase(0, 13));  // remove "extrapWeight_" and convert to int
    m_relevantLayers->push_back(layer);
  }
  return true;
}

// Streamer()
void TFCSPredictExtrapWeights::Streamer(TBuffer& buffer)
{
  // Stream an object of class TFCSPredictExtrapWeights

  if (buffer.IsReading()) {
    buffer.ReadClassBuffer(TFCSPredictExtrapWeights::Class(), this);
    m_nn.reset();
    if (m_input != nullptr && !m_input->empty()) {
      std::stringstream sin;
      sin.str(*m_input);
      auto config = lwt::parse_json(sin);
      m_nn = std::make_unique<lwt::LightweightNeuralNetwork>(
          config.inputs, config.layers, config.outputs);
    }
#ifndef __FastCaloSimStandAlone__
    // When running inside Athena, delete input/config/normInputs to free the
    // memory
    if (freemem()) {
      delete m_input;
      m_input = nullptr;
    }
#endif
  } else {
    buffer.WriteClassBuffer(TFCSPredictExtrapWeights::Class(), this);
  }
}

void TFCSPredictExtrapWeights::Print(Option_t* option) const
{
  const auto [shortprint, longprint, optprint] = print_options(option);
  TFCSLateralShapeParametrizationHitBase::Print(option);

  if (longprint) {
    FCS_MSG_INFO(optprint << "  m_input (TFCSPredictExtrapWeights): "
                          << m_input);
  }
  if (longprint) {
    FCS_MSG_INFO(optprint << "  Address of m_nn: "
                          << static_cast<const void*>(m_nn.get()));
  }
}
