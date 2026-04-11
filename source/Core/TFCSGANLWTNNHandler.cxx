// Copyright (c) 2024 CERN for the benefit of the FastCaloSim project

#include "FastCaloSim/Core/TFCSGANLWTNNHandler.h"

// For writing to a tree
#include "TBranch.h"
#include "TTree.h"

// LWTNN
#include "lwtnn/LightweightGraph.hh"
#include "lwtnn/parse_json.hh"

using namespace FastCaloSim::Core;

TFCSGANLWTNNHandler::TFCSGANLWTNNHandler(const std::string& inputFile)
    : VNetworkLWTNN(inputFile)
{
  MSG_DEBUG("Setting up from inputFile.");
  setupPersistedVariables();
  TFCSGANLWTNNHandler::setupNet();
};

TFCSGANLWTNNHandler::TFCSGANLWTNNHandler(const TFCSGANLWTNNHandler& copy_from)
    : VNetworkLWTNN(copy_from)
{
  // Cannot take copies of lwt::LightweightGraph
  // (copy constructor disabled)
  MSG_DEBUG("Making a new m_lwtnn_graph for copied network");
  std::stringstream json_stream(m_json);
  const ::lwt::GraphConfig config = ::lwt::parse_json_graph(json_stream);
  m_lwtnn_graph = std::make_unique<::lwt::LightweightGraph>(config);
  m_outputLayers = copy_from.m_outputLayers;
};

void TFCSGANLWTNNHandler::setupNet()
{
  // Backcompatability, previous versions stored this in m_input
  if (m_json.length() == 0 && m_input != nullptr) {
    m_json = *m_input;
    delete m_input;
    m_input = nullptr;
  }
  // build the graph
  MSG_VERBOSE("m_json has size " << m_json.length());
  MSG_DEBUG("m_json starts with  " << m_json.substr(0, 10));
  MSG_VERBOSE("Reading the m_json string stream into a graph network");
  std::stringstream json_stream(m_json);
  const ::lwt::GraphConfig config = ::lwt::parse_json_graph(json_stream);
  m_lwtnn_graph = std::make_unique<::lwt::LightweightGraph>(config);
  // Get the output layers
  MSG_VERBOSE("Getting output layers for neural network");
  for (auto node : config.outputs) {
    const std::string node_name = node.first;
    const ::lwt::OutputNodeConfig node_config = node.second;
    for (const std::string& label : node_config.labels) {
      MSG_VERBOSE("Found output layer called " << node_name << "_" << label);
      m_outputLayers.push_back(node_name + "_" + label);
    }
  };
  MSG_VERBOSE("Removing prefix from stored layers.");
  removePrefixes(m_outputLayers);
  MSG_VERBOSE("Finished output nodes.");
};

std::vector<std::string> TFCSGANLWTNNHandler::getOutputLayers() const
{
  return m_outputLayers;
};

// This is implement the specific compute, and ensure the output is returned in
// regular format. For LWTNN, that's easy.
TFCSGANLWTNNHandler::NetworkOutputs TFCSGANLWTNNHandler::compute(
    TFCSGANLWTNNHandler::NetworkInputs const& inputs) const
{
  MSG_DEBUG("Running computation on LWTNN graph network");
  NetworkInputs local_copy = inputs;
  if (inputs.find("Noise") != inputs.end()) {
    // Graphs from EnergyAndHitsGANV2 have the local_copy encoded as Noise =
    // node_0 and mycond = node_1
    auto noiseNode = local_copy.extract("Noise");
    noiseNode.key() = "node_0";
    local_copy.insert(std::move(noiseNode));
    auto mycondNode = local_copy.extract("mycond");
    mycondNode.key() = "node_1";
    local_copy.insert(std::move(mycondNode));
  }
  // now we can compute
  TFCSGANLWTNNHandler::NetworkOutputs outputs =
      m_lwtnn_graph->compute(local_copy);
  removePrefixes(outputs);
  MSG_DEBUG("Computation on LWTNN graph network done, returning.");
  return outputs;
};

// Giving this it's own streamer to call setupNet
void TFCSGANLWTNNHandler::Streamer(TBuffer& buf)
{
  MSG_DEBUG("In streamer of " << __FILE__);
  if (buf.IsReading()) {
    MSG_DEBUG("Reading buffer in TFCSGANLWTNNHandler ");
    // Get the persisted variables filled in
    TFCSGANLWTNNHandler::Class()->ReadBuffer(buf, this);
    MSG_DEBUG("m_json has size " << m_json.length());
    MSG_DEBUG("m_json starts with  " << m_json.substr(0, 10));
    // Setup the net, creating the non persisted variables
    // exactly as in the constructor
    this->setupNet();
#ifndef __FastCaloSimStandAlone__
    // When running inside Athena, delete persisted information
    // to conserve memory
    this->deleteAllButNet();
#endif
  } else {
    if (!m_json.empty()) {
      MSG_DEBUG("Writing buffer in TFCSGANLWTNNHandler ");
    } else {
      MSG_WARNING("Writing buffer in TFCSGANLWTNNHandler, but m_json is empty");
    };
    // Persist variables
    TFCSGANLWTNNHandler::Class()->WriteBuffer(buf, this);
  };
};
