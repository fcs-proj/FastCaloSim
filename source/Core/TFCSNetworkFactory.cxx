// Copyright (c) 2024 CERN for the benefit of the FastCaloSim project

#include <fstream>  // For checking if files exist
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
  ISF_FCS::MLogging logger;
  const std::string ending = ".*";
  const int ending_len = ending.length();
  const int filename_len = filename.length();
  if (filename_len < ending_len) {
    FCS_MSG_NOCLASS(logger, "Filename is implausably short.");
  } else if (0
             == filename.compare(filename_len - ending_len, ending_len, ending))
  {
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
};

bool TFCSNetworkFactory::isOnnxFile(std::string const& filename)
{
  ISF_FCS::MLogging logger;
  const std::string ending = ".onnx";
  const int ending_len = ending.length();
  const int filename_len = filename.length();
  bool is_onnx;
  if (filename_len < ending_len) {
    is_onnx = false;
  } else {
    is_onnx =
        (0 == filename.compare(filename_len - ending_len, ending_len, ending));
  };
  return is_onnx;
};

std::unique_ptr<VNetworkBase> TFCSNetworkFactory::create(
    std::vector<char> const& input)
{
  ISF_FCS::MLogging logger;
  FCS_MSG_NOCLASS(
      logger,
      "Directly creating ONNX network from bytes length " << input.size());
  std::unique_ptr<VNetworkBase> created(new TFCSONNXHandler(input));
  return created;
};

std::unique_ptr<VNetworkBase> TFCSNetworkFactory::create(std::string input)
{
  ISF_FCS::MLogging logger;
  resolveGlobs(input);
  if (VNetworkBase::isFile(input) && isOnnxFile(input)) {
    FCS_MSG_NOCLASS(logger,
                    "Creating ONNX network from file ..."
                        << input.substr(input.length() - 10));
    std::unique_ptr<VNetworkBase> created(new TFCSONNXHandler(input));
    return created;
  } else {
    try {
      std::unique_ptr<VNetworkBase> created(new TFCSSimpleLWTNNHandler(input));
      FCS_MSG_NOCLASS(logger,
                      "Succeeded in creating LWTNN nn from string starting "
                          << input.substr(0, 10));
      return created;
    } catch (const boost::property_tree::ptree_bad_path& e) {
      // If we get this error, it was actually a graph, not a NeuralNetwork
      std::unique_ptr<VNetworkBase> created(new TFCSGANLWTNNHandler(input));
      FCS_MSG_NOCLASS(logger, "Succeeded in creating LWTNN graph from string");
      return created;
    };
  };
};

std::unique_ptr<VNetworkBase> TFCSNetworkFactory::create(std::string input,
                                                         bool graph_form)
{
  ISF_FCS::MLogging logger;
  resolveGlobs(input);
  if (VNetworkBase::isFile(input) && isOnnxFile(input)) {
    FCS_MSG_NOCLASS(logger,
                    "Creating ONNX network from file ..."
                        << input.substr(input.length() - 10));
    std::unique_ptr<VNetworkBase> created(new TFCSONNXHandler(input));
    return created;
  } else if (graph_form) {
    FCS_MSG_NOCLASS(logger, "Creating LWTNN graph from string");
    std::unique_ptr<VNetworkBase> created(new TFCSGANLWTNNHandler(input));
    return created;
  } else {
    std::unique_ptr<VNetworkBase> created(new TFCSSimpleLWTNNHandler(input));
    FCS_MSG_NOCLASS(logger, "Creating LWTNN nn from string");
    return created;
  };
};

std::unique_ptr<VNetworkBase> TFCSNetworkFactory::create(
    std::vector<char> const& vector_input, std::string string_input)
{
  ISF_FCS::MLogging logger;
  FCS_MSG_NOCLASS(logger, "Given both bytes and a string to create an nn.");
  resolveGlobs(string_input);
  if (vector_input.size() > 0) {
    FCS_MSG_NOCLASS(logger,
                    "Bytes contains data, size=" << vector_input.size()
                                                 << ", creating from bytes.");
    return create(vector_input);
  } else if (string_input.length() > 0) {
    FCS_MSG_NOCLASS(logger,
                    "No data in bytes, string contains data, "
                        << "creating from string.");
    return create(string_input);
  } else {
    throw std::invalid_argument(
        "Neither vector_input nor string_input contain data");
  };
};

std::unique_ptr<VNetworkBase> TFCSNetworkFactory::create(
    std::vector<char> const& vector_input,
    std::string string_input,
    bool graph_form)
{
  ISF_FCS::MLogging logger;
  FCS_MSG_NOCLASS(
      logger,
      "Given both bytes, a string and graph form specified to create an nn.");
  resolveGlobs(string_input);
  if (vector_input.size() > 0) {
    FCS_MSG_NOCLASS(logger,
                    "Bytes contains data, size=" << vector_input.size()
                                                 << ", creating from bytes.");
    return create(vector_input);
  } else if (string_input.length() > 0) {
    FCS_MSG_NOCLASS(logger,
                    "No data in bytes, string contains data, "
                        << "creating from string.");
    return create(string_input, graph_form);
  } else {
    throw std::invalid_argument(
        "Neither vector_input nor string_input contain data");
  };
};
