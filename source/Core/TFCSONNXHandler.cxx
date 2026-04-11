// Copyright (c) 2024 CERN for the benefit of the FastCaloSim project

// See header file for documentation.
#include "FastCaloSim/Core/TFCSONNXHandler.h"

// For reading the binary onnx files
#include <fstream>
#include <iterator>
#include <vector>

// ONNX Runtime include(s).
#include <onnxruntime_cxx_api.h>

// For reading and writing to root
#include "TBranch.h"
#include "TFile.h"
#include "TTree.h"

// For throwing exceptions
#include <stdexcept>

using namespace FastCaloSim::Core;

TFCSONNXHandler::TFCSONNXHandler(const std::string& inputFile)
    : VNetworkBase(inputFile)
{
  MSG_INFO("Setting up from inputFile.");
  TFCSONNXHandler::setupPersistedVariables();
  TFCSONNXHandler::setupNet();
  MSG_DEBUG("Setup from file complete");
};

TFCSONNXHandler::TFCSONNXHandler(const std::vector<char>& bytes)
    : m_bytes(bytes)
{
  MSG_INFO("Given onnx session bytes as input.");
  // The super constructor got no inputFile,
  // so it won't call setupNet itself
  TFCSONNXHandler::setupNet();
  MSG_DEBUG("Setup from session complete");
};

TFCSONNXHandler::TFCSONNXHandler(const TFCSONNXHandler& copy_from)
    : VNetworkBase(copy_from)
{
  MSG_DEBUG("TFCSONNXHandler copy constructor called");
  m_bytes = copy_from.m_bytes;
  // Cannot copy a session
  // m_session = copy_from.m_session;
  // But can read it from bytes
  readSerializedSession();
  m_inputNodeNames = copy_from.m_inputNodeNames;
  m_inputNodeDims = copy_from.m_inputNodeDims;
  m_outputNodeNames = copy_from.m_outputNodeNames;
  m_outputNodeDims = copy_from.m_outputNodeDims;
  m_outputLayers = copy_from.m_outputLayers;
};

TFCSONNXHandler::NetworkOutputs TFCSONNXHandler::compute(
    TFCSONNXHandler::NetworkInputs const& inputs) const
{
  return m_computeLambda(inputs);
};

// Writing out to ttrees
void TFCSONNXHandler::writeNetToTTree(TTree& tree)
{
  MSG_DEBUG("TFCSONNXHandler writing net to tree.");
  this->writeBytesToTTree(tree, m_bytes);
};

std::vector<std::string> TFCSONNXHandler::getOutputLayers() const
{
  MSG_DEBUG("TFCSONNXHandler output layers requested.");
  return m_outputLayers;
};

void TFCSONNXHandler::deleteAllButNet()
{
  // As we don't copy the bytes, and the inputFile
  // is at most a name, nothing is needed here.
  MSG_DEBUG("Deleted nothing for ONNX.");
};

void TFCSONNXHandler::print(std::ostream& strm) const
{
  if (m_inputFile.empty()) {
    strm << "Unknown network";
  } else {
    strm << m_inputFile;
  };
  strm << "\nHas input nodes (name:dimensions);\n";
  for (size_t inp_n = 0; inp_n < m_inputNodeNames.size(); inp_n++) {
    strm << "\t" << m_inputNodeNames[inp_n] << ":[";
    for (int dim : m_inputNodeDims[inp_n]) {
      strm << " " << dim << ",";
    };
    strm << "]\n";
  };
  strm << "\nHas output nodes (name:dimensions);\n";
  for (size_t out_n = 0; out_n < m_outputNodeNames.size(); out_n++) {
    strm << "\t" << m_outputNodeNames[out_n] << ":[";
    for (int dim : m_outputNodeDims[out_n]) {
      strm << " " << dim << ",";
    };
    strm << "]\n";
  };
};

void TFCSONNXHandler::setupPersistedVariables()
{
  MSG_DEBUG("Setting up persisted variables for ONNX network.");
  // depending which constructor was called,
  // bytes may already be filled
  if (m_bytes.empty()) {
    m_bytes = getSerializedSession();
  };
  MSG_DEBUG("Setup persisted variables for ONNX network.");
};

void TFCSONNXHandler::setupNet()
{
  // From
  // https://gitlab.cern.ch/atlas/athena/-/blob/master/Control/AthOnnxruntimeUtils/AthOnnxruntimeUtils/OnnxUtils.h
  // m_session = AthONNX::CreateORTSession(inputFile);
  // This segfaults.

  // TODO; should I be using m_session_options? see
  // https://github.com/microsoft/onnxruntime-inference-examples/blob/2b42b442526b9454d1e2d08caeb403e28a71da5f/c_cxx/squeezenet/main.cpp#L71
  MSG_INFO("Setting up ONNX session.");
  this->readSerializedSession();

  // Need the type from the first node (which will be used to set
  // just set it to undefined to avoid not initialized warnings
  ONNXTensorElementDataType first_input_type =
      ONNX_TENSOR_ELEMENT_DATA_TYPE_UNDEFINED;

  // iterate over all input nodes
  MSG_DEBUG("Getting input nodes.");
  const int num_input_nodes = m_session->GetInputCount();
  Ort::AllocatorWithDefaultOptions allocator;
  for (int i = 0; i < num_input_nodes; i++) {
#if ORT_API_VERSION > 11
    Ort::AllocatedStringPtr node_names =
        m_session->GetInputNameAllocated(i, allocator);
    m_storeInputNodeNames.push_back(std::move(node_names));
    const char* input_name = m_storeInputNodeNames.back().get();
#else
    const char* input_name = m_session->GetInputName(i, allocator);
#endif
    m_inputNodeNames.push_back(input_name);
    MSG_VERBOSE("Found input node named " << input_name);

    Ort::TypeInfo type_info = m_session->GetInputTypeInfo(i);

    // For some reason unless auto is used as the return type
    // this causes a segfault once the loop ends....
    auto tensor_info = type_info.GetTensorTypeAndShapeInfo();
    if (i == 0)
      first_input_type = tensor_info.GetElementType();
    // Check the type has not changed
    if (tensor_info.GetElementType() != first_input_type) {
      MSG_ERROR("First type was " << first_input_type << ". In node " << i
                                  << " found type "
                                  << tensor_info.GetElementType());
      throw std::runtime_error(
          "Networks with varying input types not "
          "yet implemented in TFCSONNXHandler.");
    };

    std::vector<int64_t> recieved_dimension = tensor_info.GetShape();
    MSG_VERBOSE("There are " << recieved_dimension.size() << " dimensions.");
    // This vector sometimes includes a symbolic dimension
    // which is represented by -1
    // A symbolic dimension is usually a conversion error,
    // from a numpy array with a shape like (None, 7),
    // in which case it's safe to treat it as having
    // dimension 1.
    std::vector<int64_t> dimension_of_node;
    for (int64_t node_dim : recieved_dimension) {
      if (node_dim < 1) {
        MSG_WARNING("Found symbolic dimension "
                    << node_dim << " in node named " << input_name
                    << ". Will treat this as dimension 1.");
        dimension_of_node.push_back(1);
      } else {
        dimension_of_node.push_back(node_dim);
      };
    };
    m_inputNodeDims.push_back(dimension_of_node);
  };
  MSG_DEBUG("Finished looping on inputs.");

  // Outputs
  // Store the type from the first node (which will be used to set
  // m_computeLambda)
  ONNXTensorElementDataType first_output_type =
      ONNX_TENSOR_ELEMENT_DATA_TYPE_UNDEFINED;

  // iterate over all output nodes
  int num_output_nodes = m_session->GetOutputCount();
  MSG_DEBUG("Getting " << num_output_nodes << " output nodes.");
  for (int i = 0; i < num_output_nodes; i++) {
#if ORT_API_VERSION > 11
    Ort::AllocatedStringPtr node_names =
        m_session->GetOutputNameAllocated(i, allocator);
    m_storeOutputNodeNames.push_back(std::move(node_names));
    const char* output_name = m_storeOutputNodeNames.back().get();
#else
    const char* output_name = m_session->GetOutputName(i, allocator);
#endif
    m_outputNodeNames.push_back(output_name);
    MSG_VERBOSE("Found output node named " << output_name);

    const Ort::TypeInfo type_info = m_session->GetOutputTypeInfo(i);
    auto tensor_info = type_info.GetTensorTypeAndShapeInfo();
    if (i == 0)
      first_output_type = tensor_info.GetElementType();

    // Check the type has not changed
    if (tensor_info.GetElementType() != first_output_type) {
      MSG_ERROR("First type was " << first_output_type << ". In node " << i
                                  << " found type "
                                  << tensor_info.GetElementType());
      throw std::runtime_error(
          "Networks with varying output types not "
          "yet implemented in TFCSONNXHandler.");
    };

    const std::vector<int64_t> recieved_dimension = tensor_info.GetShape();
    MSG_VERBOSE("There are " << recieved_dimension.size() << " dimensions.");
    // Again, check for symbolic dimensions
    std::vector<int64_t> dimension_of_node;
    int node_size = 1;
    for (int64_t node_dim : recieved_dimension) {
      if (node_dim < 1) {
        MSG_WARNING("Found symbolic dimension "
                    << node_dim << " in node named " << output_name
                    << ". Will treat this as dimension 1.");
        dimension_of_node.push_back(1);
      } else {
        dimension_of_node.push_back(node_dim);
        node_size *= node_dim;
      };
    };
    m_outputNodeDims.push_back(dimension_of_node);
    m_outputNodeSize.push_back(node_size);

    // The outputs are treated as a flat vector
    for (int part_n = 0; part_n < node_size; part_n++) {
      // compose the output name
      std::string layer_name =
          std::string(output_name) + "_" + std::to_string(part_n);
      MSG_VERBOSE("Found output layer named " << layer_name);
      m_outputLayers.push_back(layer_name);
    }
  }
  MSG_DEBUG("Removing prefix from stored layers.");
  removePrefixes(m_outputLayers);
  MSG_DEBUG("Finished output nodes.");

  MSG_DEBUG("Setting up m_computeLambda with input type "
            << first_input_type << " and output type " << first_output_type);
  if (first_input_type == ONNX_TENSOR_ELEMENT_DATA_TYPE_FLOAT
      && first_output_type == ONNX_TENSOR_ELEMENT_DATA_TYPE_FLOAT)
  {
    // got to capture this in the lambda so it can access class methods
    m_computeLambda = [this](NetworkInputs const& inputs)
    { return computeTemplate<float, float>(inputs); };
  } else if (first_input_type == ONNX_TENSOR_ELEMENT_DATA_TYPE_DOUBLE
             && first_output_type == ONNX_TENSOR_ELEMENT_DATA_TYPE_DOUBLE)
  {
    m_computeLambda = [this](NetworkInputs const& inputs)
    { return computeTemplate<double, double>(inputs); };
  } else {
    throw std::runtime_error(
        "Haven't yet implemented that combination of "
        "input and output types as a subclass of VState.");
  };
  MSG_DEBUG("Finished setting lambda function.");
};

// Needs to also work if the input file is a root file
std::vector<char> TFCSONNXHandler::getSerializedSession(
    const std::string& tree_name)
{
  MSG_DEBUG("Getting serialized session for ONNX network.");

  if (this->isRootFile()) {
    MSG_INFO("Reading bytes from root file.");
    TFile tfile(this->m_inputFile.c_str(), "READ");
    TTree* tree = (TTree*)tfile.Get(tree_name.c_str());
    std::vector<char> bytes = this->readBytesFromTTree(*tree);
    MSG_DEBUG("Found bytes size " << bytes.size());
    return bytes;
  } else {
    MSG_INFO("Reading bytes from text file.");
    // see https://stackoverflow.com/a/50317432
    std::ifstream input(this->m_inputFile, std::ios::binary);

    std::vector<char> bytes((std::istreambuf_iterator<char>(input)),
                            (std::istreambuf_iterator<char>()));

    input.close();
    MSG_DEBUG("Found bytes size " << bytes.size());
    return bytes;
  }
};

std::vector<char> TFCSONNXHandler::readBytesFromTTree(TTree& tree)
{
  MSG_DEBUG("TFCSONNXHandler reading bytes from tree.");
  std::vector<char> bytes;
  char data;
  tree.SetBranchAddress("serialized_m_session", &data);
  for (int i = 0; tree.LoadTree(i) >= 0; i++) {
    tree.GetEntry(i);
    bytes.push_back(data);
  };
  MSG_DEBUG("TFCSONNXHandler read bytes from tree.");
  return bytes;
};

void TFCSONNXHandler::writeBytesToTTree(TTree& tree,
                                        const std::vector<char>& bytes)
{
  MSG_DEBUG("TFCSONNXHandler writing bytes to tree.");
  char m_session_data;
  tree.Branch(
      "serialized_m_session", &m_session_data, "serialized_m_session/B");
  for (Char_t here : bytes) {
    m_session_data = here;
    tree.Fill();
  };
  tree.Write();
  MSG_DEBUG("TFCSONNXHandler written bytes to tree.");
};

void TFCSONNXHandler::readSerializedSession()
{
  MSG_DEBUG("Transforming bytes to session.");
  Ort::Env env(ORT_LOGGING_LEVEL_WARNING, "test");
  Ort::SessionOptions opts;
  opts.SetInterOpNumThreads(1);
  opts.SetIntraOpNumThreads(1);
  // Prevent ONNX from spawning additional threads
  m_session =
      std::make_unique<Ort::Session>(env, m_bytes.data(), m_bytes.size(), opts);
  MSG_DEBUG("Transformed bytes to session.");
};

template<typename Tin, typename Tout>
VNetworkBase::NetworkOutputs TFCSONNXHandler::computeTemplate(
    VNetworkBase::NetworkInputs const& inputs)
{
  // working from
  // https://github.com/microsoft/onnxruntime-inference-examples/blob/main/c_cxx/squeezenet/main.cpp#L71
  //  and
  //  https://github.com/microsoft/onnxruntime-inference-examples/blob/main/c_cxx/MNIST/MNIST.cpp
  MSG_DEBUG("Setting up inputs for computation on ONNX network.");
  MSG_DEBUG("Input type " << typeid(Tin).name() << " output type "
                          << typeid(Tout).name());

  //  The inputs must be reformatted to the correct data structure.
  const size_t num_input_nodes = m_inputNodeNames.size();
  //  A pointer to all the nodes we will make
  //  Going to keep the data in each node flat, because that's easier
  std::vector<std::vector<Tin>> input_values(num_input_nodes);
  std::vector<Ort::Value> node_values;
  // Non const values that will be needed at each step.
  std::string node_name;
  int n_dimensions, elements_in_node, key_number;
  size_t first_digit;
  // Move along the list of node names gathered in the constructor
  // we need both the node name, and the dimension
  // so we cannot iterate directly on the vector.
  MSG_DEBUG("Looping over " << num_input_nodes
                            << " input nodes of ONNX network.");
  for (size_t node_n = 0; node_n < m_inputNodeNames.size(); node_n++) {
    MSG_DEBUG("Node n = " << node_n);
    node_name = m_inputNodeNames[node_n];
    MSG_DEBUG("Node name " << node_name);
    // Get the shape of this node
    n_dimensions = m_inputNodeDims[node_n].size();
    MSG_DEBUG("Node dimensions " << n_dimensions);
    elements_in_node = 1;
    for (int dimension_len : m_inputNodeDims[node_n]) {
      elements_in_node *= dimension_len;
    };
    MSG_DEBUG("Elements in node " << elements_in_node);
    for (auto inp : inputs) {
      MSG_DEBUG("Have input named " << inp.first);
    };
    // Get the node content and remove any common prefix from the elements
    const std::map<std::string, double> node_inputs = inputs.at(node_name);
    std::vector<Tin> node_elements(elements_in_node);

    MSG_DEBUG("Found node named " << node_name << " with " << elements_in_node
                                  << " elements.");
    // Then the rest should be numbers from 0 up
    for (auto element : node_inputs) {
      first_digit = element.first.find_first_of("0123456789");
      // if there is no digit, it's not an element
      if (first_digit < element.first.length()) {
        key_number = std::stoi(element.first.substr(first_digit));
        node_elements[key_number] = element.second;
      }
    }
    input_values[node_n] = node_elements;

    MSG_DEBUG("Creating ort tensor n_dimensions = "
              << n_dimensions << ", elements_in_node = " << elements_in_node);
    // Doesn't copy data internally, so vector arguments need to stay alive
    Ort::Value node =
        Ort::Value::CreateTensor<Tin>(m_memoryInfo,
                                      input_values[node_n].data(),
                                      elements_in_node,
                                      m_inputNodeDims[node_n].data(),
                                      n_dimensions);
    // Problems with the string steam when compiling separately.
    // MSG_DEBUG("Created input node " << node << " from values " <<
    // input_values[node_n]);

    node_values.push_back(std::move(node));
  }

  MSG_DEBUG("Running computation on ONNX network.");
  // All inputs have been correctly formatted and the net can be run.
  auto output_tensors = m_session->Run(Ort::RunOptions {nullptr},
                                       m_inputNodeNames.data(),
                                       &node_values[0],
                                       num_input_nodes,
                                       m_outputNodeNames.data(),
                                       m_outputNodeNames.size());

  MSG_DEBUG("Sorting outputs from computation on ONNX network.");
  // Finally, the output must be rearanged in the expected format.
  TFCSONNXHandler::NetworkOutputs outputs;
  // as the output format is just a string to double map
  // the outputs will be keyed like "<node_name>_<part_n>"
  std::string output_name;
  const Tout* output_node;
  for (size_t node_n = 0; node_n < m_outputNodeNames.size(); node_n++) {
    // get a pointer to the data
    output_node = output_tensors[node_n].GetTensorMutableData<Tout>();
    MSG_VERBOSE("output node " << output_node);
    elements_in_node = m_outputNodeSize[node_n];
    node_name = m_outputNodeNames[node_n];
    // Does the GetTensorMutableData really always return a
    // flat array?
    // Likely yes, see use of memcopy on line 301 of
    // onnxruntime/core/languge_interop_ops/pyop/pyop.cc
    for (int part_n = 0; part_n < elements_in_node; part_n++) {
      MSG_VERBOSE("Node part " << part_n << " contains "
                               << output_node[part_n]);
      // compose the output name
      output_name = node_name + "_" + std::to_string(part_n);
      outputs[output_name] = static_cast<double>(output_node[part_n]);
    }
  }
  removePrefixes(outputs);
  MSG_DEBUG("Returning outputs from computation on ONNX network.");
  return outputs;
};

// Possible to avoid copy?
// https://github.com/microsoft/onnxruntime/issues/8328
// https://github.com/microsoft/onnxruntime/pull/11789
// https://github.com/microsoft/onnxruntime/pull/8502

// Giving this its own streamer to call setupNet
void TFCSONNXHandler::Streamer(TBuffer& buf)
{
  MSG_DEBUG("In TFCSONNXHandler streamer.");
  if (buf.IsReading()) {
    MSG_INFO("Reading buffer in TFCSONNXHandler ");
    // Get the persisted variables filled in
    TFCSONNXHandler::Class()->ReadBuffer(buf, this);
    // Setup the net, creating the non persisted variables
    // exactly as in the constructor
    this->setupNet();
#ifndef __FastCaloSimStandAlone__
    // When running inside Athena, delete persisted information
    // to conserve memory
    this->deleteAllButNet();
#endif
  } else {
    MSG_INFO("Writing buffer in TFCSONNXHandler ");
    // Persist variables
    TFCSONNXHandler::Class()->WriteBuffer(buf, this);
  };
  MSG_DEBUG("Finished TFCSONNXHandler streamer.");
};
