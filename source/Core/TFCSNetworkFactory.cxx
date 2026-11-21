// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#include <fstream>  // For checking if files exist
#include <memory>
#include <stdexcept>

#include "FastCaloSim/Core/TFCSNetworkFactory.h"

#include <boost/property_tree/ptree.hpp>

#include "FastCaloSim/Core/TFCSGANLWTNNHandler.h"
#include "FastCaloSim/Core/TFCSONNXHandler.h"
#include "FastCaloSim/Core/TFCSSimpleLWTNNHandler.h"
#include "FastCaloSim/Core/VNetworkBase.h"

// For messaging
#include "FastCaloSim/Core/MLogging.h"

void TFCSNetworkFactory::resolveGlobs(std::string& filename)
{
  FastCaloSim::MLogging const logger;
  const std::string ending = ".*";
  if (filename.length() < ending.length()) {
    FCS_MSG_NOCLASS(logger, "Filename is implausably short.");
  } else if (filename.ends_with(ending)) {
    FCS_MSG_NOCLASS(logger, "Filename ends in glob.");
    // Remove the glob
    filename.pop_back();
    if (std::filesystem::exists(filename + "onnx")) {
      filename += "onnx";
    } else if (std::filesystem::exists(filename + "json")) {
      filename += std::string("json");
    } else {
      throw std::invalid_argument("No file found matching globbed filename "
                                  + filename);
    };
  };
}

auto TFCSNetworkFactory::isOnnxFile(std::string const& filename) -> bool
{
  return filename.ends_with(".onnx");
}

auto TFCSNetworkFactory::create(std::vector<char> const& input)
    -> std::unique_ptr<VNetworkBase>
{
  FastCaloSim::MLogging const logger;
  FCS_MSG_NOCLASS(
      logger,
      "Directly creating ONNX network from bytes length " << input.size());
  return std::make_unique<TFCSONNXHandler>(input);
}

auto TFCSNetworkFactory::create(std::string input)
    -> std::unique_ptr<VNetworkBase>
{
  FastCaloSim::MLogging const logger;
  resolveGlobs(input);
  if (VNetworkBase::isFile(input) && isOnnxFile(input)) {
    FCS_MSG_NOCLASS(logger,
                    "Creating ONNX network from file ..."
                        << input.substr(input.length() - 10));
    return std::make_unique<TFCSONNXHandler>(input);
  }
  try {
    auto created = std::make_unique<TFCSSimpleLWTNNHandler>(input);
    FCS_MSG_NOCLASS(logger,
                    "Succeeded in creating LWTNN nn from string starting "
                        << input.substr(0, 10));
    return created;
  } catch (const boost::property_tree::ptree_bad_path& e) {
    // If we get this error, it was actually a graph, not a NeuralNetwork
    auto created = std::make_unique<TFCSGANLWTNNHandler>(input);
    FCS_MSG_NOCLASS(logger, "Succeeded in creating LWTNN graph from string");
    return created;
  };
  ;
}

auto TFCSNetworkFactory::create(std::string input, bool graph_form)
    -> std::unique_ptr<VNetworkBase>
{
  FastCaloSim::MLogging const logger;
  resolveGlobs(input);
  if (VNetworkBase::isFile(input) && isOnnxFile(input)) {
    FCS_MSG_NOCLASS(logger,
                    "Creating ONNX network from file ..."
                        << input.substr(input.length() - 10));
    return std::make_unique<TFCSONNXHandler>(input);
  }
  if (graph_form) {
    FCS_MSG_NOCLASS(logger, "Creating LWTNN graph from string");
    return std::make_unique<TFCSGANLWTNNHandler>(input);
  } else {
    auto created = std::make_unique<TFCSSimpleLWTNNHandler>(input);
    FCS_MSG_NOCLASS(logger, "Creating LWTNN nn from string");
    return created;
  };
}

auto TFCSNetworkFactory::create(std::vector<char> const& vector_input,
                                std::string string_input)
    -> std::unique_ptr<VNetworkBase>
{
  FastCaloSim::MLogging const logger;
  FCS_MSG_NOCLASS(logger, "Given both bytes and a string to create an nn.");
  resolveGlobs(string_input);
  if (!vector_input.empty()) {
    FCS_MSG_NOCLASS(logger,
                    "Bytes contains data, size=" << vector_input.size()
                                                 << ", creating from bytes.");
    return create(vector_input);
  }
  if (string_input.length() > 0) {
    FCS_MSG_NOCLASS(
        logger,
        "No data in bytes, string contains data, " << "creating from string.");
    return create(std::move(string_input));
  } else {
    throw std::invalid_argument(
        "Neither vector_input nor string_input contain data");
  };
}

auto TFCSNetworkFactory::create(std::vector<char> const& vector_input,
                                std::string string_input,
                                bool graph_form)
    -> std::unique_ptr<VNetworkBase>
{
  FastCaloSim::MLogging const logger;
  FCS_MSG_NOCLASS(
      logger,
      "Given both bytes, a string and graph form specified to create an nn.");
  resolveGlobs(string_input);
  if (!vector_input.empty()) {
    FCS_MSG_NOCLASS(logger,
                    "Bytes contains data, size=" << vector_input.size()
                                                 << ", creating from bytes.");
    return create(vector_input);
  }
  if (string_input.length() > 0) {
    FCS_MSG_NOCLASS(
        logger,
        "No data in bytes, string contains data, " << "creating from string.");
    return create(std::move(string_input), graph_form);
  } else {
    throw std::invalid_argument(
        "Neither vector_input nor string_input contain data");
  };
}
