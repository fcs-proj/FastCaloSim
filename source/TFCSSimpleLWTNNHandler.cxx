/*
  Copyright (C) 2002-2024 CERN for the benefit of the ATLAS collaboration
*/

#include "FastCaloSim/TFCSSimpleLWTNNHandler.h"

// For writing to a tree
#include "TBranch.h"
#include "TTree.h"

// LWTNN
#include "lwtnn/LightweightNeuralNetwork.hh"
#include "lwtnn/parse_json.hh"

TFCSSimpleLWTNNHandler::TFCSSimpleLWTNNHandler(const std::string& inputFile)
    : VNetworkLWTNN(inputFile)
{
  ATH_MSG_DEBUG("Setting up from inputFile.");
  setupPersistedVariables();
  TFCSSimpleLWTNNHandler::setupNet();
};

TFCSSimpleLWTNNHandler::TFCSSimpleLWTNNHandler(
    const TFCSSimpleLWTNNHandler& copy_from)
    : VNetworkLWTNN(copy_from)
{
  // Cannot take copy of lwt::LightweightNeuralNetwork
  // (copy constructor disabled)
  ATH_MSG_DEBUG("Making new m_lwtnn_neural for copy of network.");
  std::stringstream json_stream(m_json);
  const lwt::JSONConfig config = lwt::parse_json(json_stream);
  m_lwtnn_neural = std::make_unique<lwt::LightweightNeuralNetwork>(
      config.inputs, config.layers, config.outputs);
  m_outputLayers = copy_from.m_outputLayers;
};

void TFCSSimpleLWTNNHandler::setupNet()
{
  // build the graph
  ATH_MSG_DEBUG("Reading the m_json string stream into a neural network");
  std::stringstream json_stream(m_json);
  const lwt::JSONConfig config = lwt::parse_json(json_stream);
  m_lwtnn_neural = std::make_unique<lwt::LightweightNeuralNetwork>(
      config.inputs, config.layers, config.outputs);
  // Get the output layers
  ATH_MSG_DEBUG("Getting output layers for neural network");
  for (std::string name : config.outputs) {
    ATH_MSG_VERBOSE("Found output layer called " << name);
    m_outputLayers.push_back(name);
  };
  ATH_MSG_DEBUG("Removing prefix from stored layers.");
  removePrefixes(m_outputLayers);
  ATH_MSG_DEBUG("Finished output nodes.");
}

std::vector<std::string> TFCSSimpleLWTNNHandler::getOutputLayers() const
{
  return m_outputLayers;
};

// This is implement the specific compute, and ensure the output is returned in
// regular format. For LWTNN, that's easy.
TFCSSimpleLWTNNHandler::NetworkOutputs TFCSSimpleLWTNNHandler::compute(
    TFCSSimpleLWTNNHandler::NetworkInputs const& inputs) const
{
  ATH_MSG_DEBUG("Running computation on LWTNN neural network");
  ATH_MSG_DEBUG(VNetworkBase::representNetworkInputs(inputs, 20));
  // Flatten the map depth
  if (inputs.size() != 1) {
    ATH_MSG_ERROR("The inputs have multiple elements."
                  << " An LWTNN neural network can only handle one node.");
  };
  std::map<std::string, double> flat_inputs;
  for (auto node : inputs) {
    flat_inputs = node.second;
  }
  // Now we have flattened, we can compute.
  NetworkOutputs outputs = m_lwtnn_neural->compute(flat_inputs);
  removePrefixes(outputs);
  ATH_MSG_DEBUG(VNetworkBase::representNetworkOutputs(outputs, 20));
  ATH_MSG_DEBUG("Computation on LWTNN neural network done, returning");
  return outputs;
};

// Giving this it's own streamer to call setupNet
void TFCSSimpleLWTNNHandler::Streamer(TBuffer& buf)
{
  ATH_MSG_DEBUG("In streamer of " << __FILE__);
  if (buf.IsReading()) {
    ATH_MSG_DEBUG("Reading buffer in TFCSSimpleLWTNNHandler ");
    // Get the persisted variables filled in
    TFCSSimpleLWTNNHandler::Class()->ReadBuffer(buf, this);
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
      ATH_MSG_DEBUG("Writing buffer in TFCSSimpleLWTNNHandler ");
    } else {
      ATH_MSG_WARNING(
          "Writing buffer in TFCSSimpleLWTNNHandler, but m_json is empty.");
    }
    // Persist variables
    TFCSSimpleLWTNNHandler::Class()->WriteBuffer(buf, this);
  };
};
