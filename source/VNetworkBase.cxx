/*
  Copyright (C) 2002-2024 CERN for the benefit of the ATLAS collaboration
*/

#include <iostream>

#include "FastCaloSim/VNetworkBase.h"

// For streamer
#include "TBuffer.h"

// For reading and writing to root
#include "TFile.h"
#include "TTree.h"

// Probably called by a streamer.
VNetworkBase::VNetworkBase()
    : m_inputFile("unknown") {};

// record the input file and provided it's not empty call SetUp
VNetworkBase::VNetworkBase(const std::string& inputFile)
    : m_inputFile(inputFile)
{
  ATH_MSG_DEBUG("Constructor called with inputFile");
};

// No setupPersistedVariables or setupNet here!
VNetworkBase::VNetworkBase(const VNetworkBase& copy_from)
    : MLogging()
    , m_inputFile(copy_from.m_inputFile) {};

// Nothing is needed from the destructor right now.
// We don't use new anywhere, so the whole thing should clean
// itself up.
VNetworkBase::~VNetworkBase() {};

std::string VNetworkBase::representNetworkInputs(
    VNetworkBase::NetworkInputs const& inputs, int maxValues)
{
  std::string representation =
      "NetworkInputs, outer size " + std::to_string(inputs.size());
  int valuesIncluded = 0;
  for (const auto& outer : inputs) {
    representation += "\n key->" + outer.first + "; ";
    for (const auto& inner : outer.second) {
      representation += inner.first + "=" + std::to_string(inner.second) + ", ";
      ++valuesIncluded;
      if (valuesIncluded > maxValues)
        break;
    };
    if (valuesIncluded > maxValues)
      break;
  };
  representation += "\n";
  return representation;
};

std::string VNetworkBase::representNetworkOutputs(
    VNetworkBase::NetworkOutputs const& outputs, int maxValues)
{
  std::string representation =
      "NetworkOutputs, size " + std::to_string(outputs.size()) + "; \n";
  int valuesIncluded = 0;
  for (const auto& item : outputs) {
    representation += item.first + "=" + std::to_string(item.second) + ", ";
    ++valuesIncluded;
    if (valuesIncluded > maxValues)
      break;
  };
  representation += "\n";
  return representation;
};

// this is also used for the stream operator
void VNetworkBase::print(std::ostream& strm) const
{
  if (m_inputFile.empty()) {
    ATH_MSG_DEBUG("Making a network without a named inputFile");
    strm << "Unknown network";
  } else {
    ATH_MSG_DEBUG("Making a network with input file " << m_inputFile);
    strm << m_inputFile;
  };
};

void VNetworkBase::writeNetToTTree(TFile& root_file,
                                   std::string const& tree_name)
{
  ATH_MSG_DEBUG("Making tree name " << tree_name);
  root_file.cd();
  const std::string title = "onnxruntime saved network";
  TTree tree(tree_name.c_str(), title.c_str());
  this->writeNetToTTree(tree);
  root_file.Write();
};

void VNetworkBase::writeNetToTTree(std::string const& root_name,
                                   std::string const& tree_name)
{
  ATH_MSG_DEBUG("Making or updating file name " << root_name);
  TFile root_file(root_name.c_str(), "UPDATE");
  this->writeNetToTTree(root_file, tree_name);
  root_file.Close();
};

bool VNetworkBase::isRootFile(std::string const& filename) const
{
  const std::string* to_check = &filename;
  if (filename.length() == 0) {
    to_check = &this->m_inputFile;
    ATH_MSG_DEBUG("No file name given, so using m_inputFile, " << m_inputFile);
  };
  const std::string ending = ".root";
  const int ending_len = ending.length();
  const int filename_len = to_check->length();
  if (filename_len < ending_len) {
    return false;
  }
  return (0
          == to_check->compare(filename_len - ending_len, ending_len, ending));
};

bool VNetworkBase::isFile() const
{
  return isFile(m_inputFile);
};

bool VNetworkBase::isFile(std::string const& inputFile)
{
  if (FILE* file = std::fopen(inputFile.c_str(), "r")) {
    std::fclose(file);
    return true;
  } else {
    return false;
  };
};

namespace
{
int GetPrefixLength(const std::vector<std::string>& strings)
{
  const std::string first = strings[0];
  int length = first.length();
  for (const std::string& this_string : strings) {
    for (int i = 0; i < length; i++) {
      if (first[i] != this_string[i]) {
        length = i;
        break;
      }
    }
  }
  return length;
};
}  // namespace

void VNetworkBase::removePrefixes(std::vector<std::string>& output_names) const
{
  const int length = GetPrefixLength(output_names);
  for (long unsigned int i = 0; i < output_names.size(); i++)
    output_names[i] = output_names[i].substr(length);
};

void VNetworkBase::removePrefixes(VNetworkBase::NetworkOutputs& outputs) const
{
  std::vector<std::string> output_layers;
  for (auto const& output : outputs)
    output_layers.push_back(output.first);
  const int length = GetPrefixLength(output_layers);
  for (std::string layer_name : output_layers) {
    // remove this output
    auto nodeHandle = outputs.extract(layer_name);
    // change the key
    nodeHandle.key() = layer_name.substr(length);
    // replace the output
    outputs.insert(std::move(nodeHandle));
  }
};
