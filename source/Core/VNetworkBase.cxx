// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#include <algorithm>
#include <cstddef>
#include <format>
#include <iostream>

#include "FastCaloSim/Core/VNetworkBase.h"

// For streamer
#include "TBuffer.h"

// For reading and writing to root
#include "TFile.h"
#include "TTree.h"

// Probably called by a streamer.
VNetworkBase::VNetworkBase()
    : m_inputFile("unknown")
{
}

// record the input file and provided it's not empty call SetUp
// NOLINTNEXTLINE(modernize-pass-by-value): keep the public signature
VNetworkBase::VNetworkBase(const std::string& inputFile)
    : m_inputFile(inputFile)
{
  FCS_MSG_DEBUG("Constructor called with inputFile");
}

// No setupPersistedVariables or setupNet here! The logging settings are not
// copied either.
// NOLINTNEXTLINE(bugprone-copy-constructor-init): default logging settings
VNetworkBase::VNetworkBase(const VNetworkBase& copy_from)
    : m_inputFile(copy_from.m_inputFile)
{
}

// Nothing is needed from the destructor right now.
// We don't use new anywhere, so the whole thing should clean
// itself up.
VNetworkBase::~VNetworkBase() = default;

auto VNetworkBase::representNetworkInputs(
    VNetworkBase::NetworkInputs const& inputs, int maxValues) -> std::string
{
  std::string representation =
      "NetworkInputs, outer size " + std::to_string(inputs.size());
  int valuesIncluded = 0;
  for (const auto& outer : inputs) {
    representation += "\n key->" + outer.first + "; ";
    for (const auto& inner : outer.second) {
      representation += std::format("{}={:f}, ", inner.first, inner.second);
      ++valuesIncluded;
      if (valuesIncluded > maxValues) {
        break;
      }
    };
    if (valuesIncluded > maxValues) {
      break;
    }
  };
  representation += "\n";
  return representation;
}

auto VNetworkBase::representNetworkOutputs(
    VNetworkBase::NetworkOutputs const& outputs, int maxValues) -> std::string
{
  std::string representation =
      "NetworkOutputs, size " + std::to_string(outputs.size()) + "; \n";
  int valuesIncluded = 0;
  for (const auto& item : outputs) {
    representation += std::format("{}={:f}, ", item.first, item.second);
    ++valuesIncluded;
    if (valuesIncluded > maxValues) {
      break;
    }
  };
  representation += "\n";
  return representation;
}

// this is also used for the stream operator
void VNetworkBase::print(std::ostream& strm) const
{
  if (m_inputFile.empty()) {
    FCS_MSG_DEBUG("Making a network without a named inputFile");
    strm << "Unknown network";
  } else {
    FCS_MSG_DEBUG("Making a network with input file " << m_inputFile);
    strm << m_inputFile;
  };
}

void VNetworkBase::writeNetToTTree(TFile& root_file,
                                   std::string const& tree_name)
{
  FCS_MSG_DEBUG("Making tree name " << tree_name);
  root_file.cd();
  const std::string title = "onnxruntime saved network";
  TTree tree(tree_name.c_str(), title.c_str());
  this->writeNetToTTree(tree);
  root_file.Write();
}

void VNetworkBase::writeNetToTTree(std::string const& root_name,
                                   std::string const& tree_name)
{
  FCS_MSG_DEBUG("Making or updating file name " << root_name);
  TFile root_file(root_name.c_str(), "UPDATE");
  this->writeNetToTTree(root_file, tree_name);
  root_file.Close();
}

auto VNetworkBase::isRootFile(std::string const& filename) const -> bool
{
  const std::string* to_check = &filename;
  if (filename.empty()) {
    to_check = &this->m_inputFile;
    FCS_MSG_DEBUG("No file name given, so using m_inputFile, " << m_inputFile);
  };
  return to_check->ends_with(".root");
}

auto VNetworkBase::isFile() const -> bool
{
  return isFile(m_inputFile);
}

auto VNetworkBase::isFile(std::string const& inputFile) -> bool
{
  if (FILE* file = std::fopen(inputFile.c_str(), "r")) {
    static_cast<void>(std::fclose(file));
    return true;
  }
  return false;
  ;
}

namespace
{
auto GetPrefixLength(const std::vector<std::string>& strings) -> std::size_t
{
  if (strings.empty()) {
    return 0;
  }
  const std::string& first = strings[0];
  std::size_t length = first.length();
  for (const std::string& this_string : strings) {
    const std::size_t max_length = std::min(length, this_string.length());
    std::size_t i = 0;
    while (i < max_length && first[i] == this_string[i]) {
      ++i;
    }
    length = i;
  }
  return length;
}
}  // namespace

void VNetworkBase::removePrefixes(std::vector<std::string>& output_names)
{
  const std::size_t length = GetPrefixLength(output_names);
  for (auto& output_name : output_names) {
    output_name = output_name.substr(length);
  }
}

void VNetworkBase::removePrefixes(VNetworkBase::NetworkOutputs& outputs)
{
  std::vector<std::string> output_layers;
  for (auto const& output : outputs) {
    output_layers.push_back(output.first);
  }
  const std::size_t length = GetPrefixLength(output_layers);
  for (const std::string& layer_name : output_layers) {
    // remove this output
    auto nodeHandle = outputs.extract(layer_name);
    // change the key
    nodeHandle.key() = layer_name.substr(length);
    // replace the output
    outputs.insert(std::move(nodeHandle));
  }
}
