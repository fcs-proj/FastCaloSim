/*
  Copyright (C) 2002-2024 CERN for the benefit of the ATLAS collaboration
*/

#include <fstream>
#include <iostream>

#include "FastCaloSim/Core/TFCSPredictExtrapWeights.h"

#include <CLHEP/Random/RanluxEngine.h>

#include "FastCaloSim/Core/TFCSExtrapolationState.h"
#include "FastCaloSim/Core/TFCSSimulationState.h"
#include "FastCaloSim/Core/TFCSTruthState.h"
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
  if (m_input != nullptr) {
    delete m_input;
  }
  if (m_relevantLayers != nullptr) {
    delete m_relevantLayers;
  }
  if (m_normLayers != nullptr) {
    delete m_normLayers;
  }
  if (m_normMeans != nullptr) {
    delete m_normMeans;
  }
  if (m_normStdDevs != nullptr) {
    delete m_normStdDevs;
  }
  if (m_nn != nullptr) {
    delete m_nn;
  }
}

bool TFCSPredictExtrapWeights::operator==(
    const TFCSParametrizationBase& ref) const
{
  if (IsA() != ref.IsA()) {
    ATH_MSG_DEBUG("operator==: different class types "
                  << IsA()->GetName() << " != " << ref.IsA()->GetName());
    return false;
  }
  const TFCSPredictExtrapWeights& ref_typed =
      static_cast<const TFCSPredictExtrapWeights&>(ref);

  if (TFCSParametrizationBase::compare(ref))
    return true;
  if (!TFCSParametrization::compare(ref))
    return false;
  if (!TFCSLateralShapeParametrization::compare(ref))
    return false;

  return (m_input->compare(*ref_typed.m_input) == 0);
}

// getNormInputs()
// Get values needed to normalize inputs
bool TFCSPredictExtrapWeights::getNormInputs(
    const std::string& etaBin, const std::string& FastCaloTXTInputFolderName)
{
  ATH_MSG_DEBUG(" Getting normalization inputs... ");

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
  std::string inputFileName = FastCaloTXTInputFolderName
      + "MeanStdDevEnergyFractions_eta_" + etaBin + ".txt";
  ATH_MSG_DEBUG(" Opening " << inputFileName);
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
            int index = std::stoi(substr.substr(substr.find('_') + 1));
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
    ATH_MSG_ERROR(" Unable to open file " << inputFileName);
    return false;
  }

  return true;
}

// prepareInputs()
// Prepare input variables to the Neural Network
std::map<std::string, double> TFCSPredictExtrapWeights::prepareInputs(
    TFCSSimulationState& simulstate, const float truthE) const
{
  std::map<std::string, double> inputVariables;
  for (int ilayer = 0; ilayer < m_geo->n_layers(); ++ilayer) {
    if (std::find(m_relevantLayers->cbegin(), m_relevantLayers->cend(), ilayer)
        != m_relevantLayers->cend())
    {
      const std::string layer = std::to_string(ilayer);
      // Find index
      auto itr =
          std::find(m_normLayers->cbegin(), m_normLayers->cend(), ilayer);
      if (itr != m_normLayers->cend()) {
        const int index = std::distance(m_normLayers->cbegin(), itr);
        inputVariables["ef_" + layer] =
            (simulstate.Efrac(ilayer) - std::as_const(m_normMeans)->at(index))
            / std::as_const(m_normStdDevs)->at(index);
      } else {
        ATH_MSG_ERROR("Normalization information not found for layer "
                      << ilayer);
      }
    }
  }
  // Find index for truth energy
  auto itr = std::find(m_normLayers->cbegin(), m_normLayers->cend(), -1);
  int index = std::distance(m_normLayers->cbegin(), itr);
  inputVariables["etrue"] = (truthE - std::as_const(m_normMeans)->at(index))
      / std::as_const(m_normStdDevs)->at(index);
  if (is_match_pdgid(22)) {
    inputVariables["pdgId"] = 1;  // one hot encoding
  } else if (is_match_pdgid(11) || is_match_pdgid(-11)) {
    inputVariables["pdgId"] = 0;  // one hot encoding
  }
  return inputVariables;
}

// simulate()
// get predicted extrapolation weights and save them as AuxInfo in simulstate
FCSReturnCode TFCSPredictExtrapWeights::simulate(
    TFCSSimulationState& simulstate,
    const TFCSTruthState* truth,
    const TFCSExtrapolationState* /*extrapol*/) const
{
  // Get inputs to Neural Network
  std::map<std::string, double> inputVariables =
      prepareInputs(simulstate, truth->E() * 0.001);

  // Get predicted extrapolation weights
  auto outputs = m_nn->compute(inputVariables);
  for (int ilayer = 0; ilayer < m_geo->n_layers(); ++ilayer) {
    if (std::find(m_relevantLayers->cbegin(), m_relevantLayers->cend(), ilayer)
        != m_relevantLayers->cend())
    {
      ATH_MSG_DEBUG("TFCSPredictExtrapWeights::simulate: layer: "
                    << ilayer << " weight: "
                    << outputs["extrapWeight_" + std::to_string(ilayer)]);
      float weight = outputs["extrapWeight_" + std::to_string(ilayer)];
      // Protections
      if (weight < 0) {
        weight = 0;
      } else if (weight > 1) {
        weight = 1;
      }
      simulstate.setAuxInfo<float>(ilayer, weight);
    } else {  // use weight=0.5 for non-relevant layers
      ATH_MSG_DEBUG(
          "Setting weight=0.5 for layer = " << std::to_string(ilayer));
      simulstate.setAuxInfo<float>(ilayer, float(0.5));
    }
  }
  return FCSSuccess;
}

// simulate_hit()
FCSReturnCode TFCSPredictExtrapWeights::simulate_hit(
    Hit& hit,
    TFCSSimulationState& simulstate,
    const TFCSTruthState* /*truth*/,
    const TFCSExtrapolationState* extrapol)
{
  const int cs = calosample();

  // Get corresponding predicted extrapolation weight from simulstate
  float extrapWeight;
  if (simulstate.hasAuxInfo(cs)) {
    extrapWeight = simulstate.getAuxInfo<float>(cs);
  } else {  // missing AuxInfo
    ATH_MSG_FATAL(
        "Simulstate is not decorated with extrapolation weights for cs = "
        << std::to_string(cs));
    return FCSFatal;
  }

  double eta = (1. - extrapWeight) * extrapol->eta(cs, Cell::SubPos::ENT)
      + extrapWeight * extrapol->eta(cs, Cell::SubPos::EXT);
  double phi = (1. - extrapWeight) * extrapol->phi(cs, Cell::SubPos::ENT)
      + extrapWeight * extrapol->phi(cs, Cell::SubPos::EXT);
  float extrapWeight_for_r_z = extrapWeight;
  if (UseHardcodedWeight()) {
    extrapWeight_for_r_z = 0.5;
    ATH_MSG_DEBUG(
        "Will use extrapWeight=0.5 for r and z when constructing a hit");
  } else {
    ATH_MSG_DEBUG(
        "Will use predicted extrapWeight also for r and z when "
        "constructing a hit");
  }
  double r = (1. - extrapWeight_for_r_z) * extrapol->r(cs, Cell::SubPos::ENT)
      + extrapWeight_for_r_z * extrapol->r(cs, Cell::SubPos::EXT);
  double z = (1. - extrapWeight_for_r_z) * extrapol->z(cs, Cell::SubPos::ENT)
      + extrapWeight_for_r_z * extrapol->z(cs, Cell::SubPos::EXT);

  if (!std::isfinite(r) || !std::isfinite(z) || !std::isfinite(eta)
      || !std::isfinite(phi))
  {
    ATH_MSG_WARNING(
        "Extrapolator contains NaN or infinite number.\nSetting "
        "center position to calo boundary.");
    ATH_MSG_WARNING("Before fix: center_r: "
                    << r << " center_z: " << z << " center_phi: " << phi
                    << " center_eta: " << eta << " weight: " << extrapWeight
                    << " cs: " << cs);
    // If extrapolator fails we can set position to calo boundary
    r = extrapol->IDCaloBoundary_r();
    z = extrapol->IDCaloBoundary_z();
    eta = extrapol->IDCaloBoundary_eta();
    phi = extrapol->IDCaloBoundary_phi();
    ATH_MSG_WARNING("After fix: center_r: "
                    << r << " center_z: " << z << " center_phi: " << phi
                    << " center_eta: " << eta << " weight: " << extrapWeight
                    << " cs: " << cs);
  }

  hit.setCenter_r(r);
  hit.setCenter_z(z);
  hit.setCenter_eta(eta);
  hit.setCenter_phi(phi);

  ATH_MSG_DEBUG("TFCSPredictExtrapWeights: center_r: "
                << hit.center_r() << " center_z: " << hit.center_z()
                << " center_phi: " << hit.center_phi()
                << " center_eta: " << hit.center_eta()
                << " weight: " << extrapWeight << " cs: " << cs);

  return FCSSuccess;
}

// initializeNetwork()
// Initialize lwtnn network
bool TFCSPredictExtrapWeights::initializeNetwork(
    int pid,
    const std::string& etaBin,
    const std::string& FastCaloNNInputFolderName)
{
  ATH_MSG_INFO(
      "Using FastCaloNNInputFolderName: " << FastCaloNNInputFolderName);
  set_pdgid(pid);

  std::string inputFileName =
      FastCaloNNInputFolderName + "NN_" + etaBin + ".json";
  ATH_MSG_DEBUG("Will read JSON file: " << inputFileName);
  if (inputFileName.empty()) {
    ATH_MSG_ERROR("Could not find json file " << inputFileName);
    return false;
  } else {
    ATH_MSG_INFO("For pid: " << pid << " and etaBin" << etaBin
                             << ", loading json file " << inputFileName);
    std::ifstream input(inputFileName);
    std::stringstream sin;
    sin << input.rdbuf();
    input.close();
    auto config = lwt::parse_json(sin);
    m_nn = new lwt::LightweightNeuralNetwork(
        config.inputs, config.layers, config.outputs);
    if (m_nn == nullptr) {
      ATH_MSG_ERROR("Could not create LightWeightNeuralNetwork from "
                    << inputFileName);
      return false;
    }
    if (m_input != nullptr) {
      delete m_input;
    }
    m_input = new std::string(sin.str());
    // Extract relevant layers from the outputs
    m_relevantLayers = new std::vector<int>();
    for (auto name : config.outputs) {
      int layer = std::stoi(
          name.erase(0, 13));  // remove "extrapWeight_" and convert to int
      m_relevantLayers->push_back(layer);
    }
  }
  return true;
}

// Streamer()
void TFCSPredictExtrapWeights::Streamer(TBuffer& R__b)
{
  // Stream an object of class TFCSPredictExtrapWeights

  if (R__b.IsReading()) {
    R__b.ReadClassBuffer(TFCSPredictExtrapWeights::Class(), this);
    if (m_nn != nullptr) {
      delete m_nn;
      m_nn = nullptr;
    }
    if (m_input && !m_input->empty()) {
      std::stringstream sin;
      sin.str(*m_input);
      auto config = lwt::parse_json(sin);
      m_nn = new lwt::LightweightNeuralNetwork(
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
    R__b.WriteClassBuffer(TFCSPredictExtrapWeights::Class(), this);
  }
}

void TFCSPredictExtrapWeights::Print(Option_t* option) const
{
  TString opt(option);
  bool shortprint = opt.Index("short") >= 0;
  bool longprint = msgLvl(MSG::DEBUG) || (msgLvl(MSG::INFO) && !shortprint);
  TString optprint = opt;
  optprint.ReplaceAll("short", "");
  TFCSLateralShapeParametrizationHitBase::Print(option);

  if (longprint)
    ATH_MSG_INFO(optprint << "  m_input (TFCSPredictExtrapWeights): "
                          << m_input);
  if (longprint)
    ATH_MSG_INFO(optprint << "  Address of m_nn: " << (void*)m_nn);
}
