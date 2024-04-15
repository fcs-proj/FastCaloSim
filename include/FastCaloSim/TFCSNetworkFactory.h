/**
 * Copyright (C) 2002-2024 CERN for the benefit of the ATLAS collaboration
 *
 * Class to perform runtime selection from the derived
 * classes of VNetworkBase given input for a network.
 *
 * Has only static functions because no statelike
 * information is needed to make this decision.
 *
 * Information about the which network would be
 * apropreate can be specified, or left entirely
 * to the factory to determine.
 *
 **/
#ifndef TFCSNETWORKFACTORY_H
#define TFCSNETWORKFACTORY_H

#include <filesystem>
#include <memory>
#include <string>
#include <vector>

#include "FastCaloSim/VNetworkBase.h"

class TFCSNetworkFactory
{
public:
  // Unspecified TFCSGANLWTNNHandler or TFCSSimpleLWTNNHandler, take a guess
  /**
   * @brief Given a string, make a network.
   *
   * This function will first check if the string is the path of a readable
   * file on the disk. If so, the file suffix is used to decide if it's an
   * onnx (.onnx) or lwtnn file (.json). If the filepath ends in .* then
   * first an onnx then an lwtnn file will be tried. The file is read and
   * parsed into a network.
   *
   * If the string is not a filepath, it is assumed to be the content of a
   * json file to make an lwtnn network.
   *
   * When an lwtnn network is created, first the TFCSSimpleLWTNNHandler
   * format is tried, and if this raises an error, the TFCSGANLWTNNHandler
   * is applied. The former is simpler than the latter, so it will always
   * fail to parse the more complex graph format.
   *
   * @param input  Either a file path, or the content of a file.
   *
   **/
  static std::unique_ptr<VNetworkBase> create(std::string input);
  // Specified TFCSGANLWTNNHandler or TFCSSimpleLWTNNHandler
  /**
   * @brief Given a string, and information about format, make a network.
   *
   * This function will first check if the string is the path of a readable
   * file on the disk. If so, the file suffix is used to decide if it's an
   * onnx (.onnx) or lwtnn file (.json). If the filepath ends in .* then
   * first an onnx then an lwtnn file will be tried. The file is read and
   * parsed into a network.
   *
   * If the string is not a filepath, it is assumed to be the content of a
   * json file to make an lwtnn network.
   *
   * When an lwtnn network is created, if graph_form is true
   * the network will be a TFCSSimpleLWTNNHandler otherwise it is
   * a TFCSGANLWTNNHandler.
   *
   * @param input        Either a file path, or the content of a file.
   * @param graph_form   Is the network the more complex graph form?
   *
   **/
  static std::unique_ptr<VNetworkBase> create(std::string input,
                                              bool graph_form);

  /**
   * @brief Given a vector of chars (bytes), make a network.
   *
   * This function will always create a TFCSONNXHandler.
   * Caution: this function is designed to modify its input.
   *
   * @param input        The content of an onnx proto file.
   *
   **/
  static std::unique_ptr<VNetworkBase> create(std::vector<char> const& input);

  /**
   * @brief Create a network from whichever input isn't empty.
   *
   * If the vector_input is not empty, construct a network from that,
   * otherwise, use the string_input to construct a network.
   *
   * @param vector_input        The content of an onnx proto file.
   * @param string_input        Either a file path, or the content of a file.
   **/
  static std::unique_ptr<VNetworkBase> create(
      std::vector<char> const& vector_input, std::string string_input);
  /**
   * @brief Create a network from whichever input isn't empty.
   *
   * If the vector_input is not empty, construct a network from that,
   * otherwise, use the string_input to construct a network.
   * Whether the network is in graph form is specified for LWTNN networks.
   *
   * @param vector_input        The content of an onnx proto file.
   * @param string_input        Either a file path, or the content of a file.
   * @param graph_form   Is the network the more complex graph form?
   **/
  static std::unique_ptr<VNetworkBase> create(
      std::vector<char> const& vector_input,
      std::string string_input,
      bool graph_form);

private:
  /**
   * @brief If the filepath ends in .* change it to .onnx or .json
   *
   * If the filepath doesn't end in .*, no change is made.
   * Will check first for a .onnx file, then look for a .json.
   * Throws an exception if neither are found.
   *
   * @param filename  Path to check.
   **/
  static void resolveGlobs(std::string& filename);

  /**
   * @brief Check if a filename seems to be an onnx file.
   *
   * Really just checks if the input ends in ".onnx"
   *
   * @param filename  Path to check.
   **/
  static bool isOnnxFile(std::string const& filename);
};

#endif  // TFCSNETWORKFACTORY_H
