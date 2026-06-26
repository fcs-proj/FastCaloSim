// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#include "FastCaloSim/Core/TFCSMLCalorimeterSimulator.h"

#include "CLHEP/Random/RandGauss.h"
#include "FastCaloSim/Core/TFCSNetworkFactory.h"

TFCSMLCalorimeterSimulator::TFCSMLCalorimeterSimulator() {}

TFCSMLCalorimeterSimulator::~TFCSMLCalorimeterSimulator() {}

bool TFCSMLCalorimeterSimulator::loadSimulator(std::string filename)
{
  // Load the simulator
  try {
    m_onnx_model = TFCSNetworkFactory::create(filename);
  } catch (std::exception& e) {
    FCS_MSG_ERROR("Failed to load simulator from file "
                  << filename << " with error " << e.what());
    return false;
  }

  if (m_onnx_model == nullptr) {
    FCS_MSG_ERROR("Failed to load simulator from file " << filename);
    return false;
  }

  return true;
}

VNetworkBase::NetworkOutputs TFCSMLCalorimeterSimulator::predictVoxels(
    TFCSSimulationState& simulstate, float eta, float energy) const
{
  // Bring eta into the needed range
  eta = std::abs(eta) * 10;

  // Initialize the energy and eta input vectors
  std::vector<float> eta_vector(m_nEvents, eta);
  std::vector<float> energy_vector(m_nEvents, energy);

  // sample the z vectors according to a standard normal distribution
  std::vector<float> z_shape_vector(m_nEvents * m_nVoxels, 0.0);
  std::vector<float> z_energy_vector(m_nEvents * m_nLayers, 0.0);
  for (auto& z_shape : z_shape_vector) {
    z_shape = CLHEP::RandGauss::shoot(simulstate.randomEngine(), 0.0, 1.0);
  }
  for (auto& z_energy : z_energy_vector) {
    z_energy = CLHEP::RandGauss::shoot(simulstate.randomEngine(), 0.0, 1.0);
  }

  // Prepare the inputs for the network
  VNetworkBase::NetworkInputs inputs;

  int i = 0;
  for (float eta : eta_vector) {
    inputs["inn_eta_in"].insert(
        std::pair<std::string, double>("variable_" + std::to_string(i), eta));
    i++;
  }

  i = 0;
  for (float energy : energy_vector) {
    inputs["inn_einc_in"].insert(std::pair<std::string, double>(
        "variable_" + std::to_string(i), energy));
    i++;
  }

  i = 0;
  for (float z_shape : z_shape_vector) {
    inputs["cfm_z_shape"].insert(std::pair<std::string, double>(
        "variable_" + std::to_string(i), z_shape));
    i++;
  }

  i = 0;
  for (float z_energy : z_energy_vector) {
    inputs["inn_z_energy"].insert(std::pair<std::string, double>(
        "variable_" + std::to_string(i), z_energy));
    i++;
  }

  // Compute the network outputs
  if (!m_onnx_model) {
    FCS_MSG_ERROR("predictVoxels called before a simulator was loaded");
    return {};
  }
  VNetworkBase::NetworkOutputs outputs = m_onnx_model->compute(inputs);

  return outputs;
}

TFCSMLCalorimeterSimulator::event_t TFCSMLCalorimeterSimulator::getEvent(
    TFCSSimulationState& simulstate, float eta, float energy) const
{
  // Get the voxel energies
  VNetworkBase::NetworkOutputs outputs = predictVoxels(simulstate, eta, energy);

  if (outputs.empty()) {
    FCS_MSG_ERROR(
        "predictVoxels returned empty outputs; returning empty event");
    return event_t {};
  }

  // Check if any output contains NaN; retry up to 5 times if so
  auto has_nan = [](const VNetworkBase::NetworkOutputs& out)
  {
    for (const auto& kv : out) {
      if (std::isnan(kv.second))
        return true;
    }
    return false;
  };

  int retry = 0;
  while (has_nan(outputs)) {
    if (retry >= 5) {
      FCS_MSG_WARNING("Network output still contains NaN after "
                      << retry << " retries. Giving up.");
      break;
    }
    FCS_MSG_WARNING("Network output contains NaN. Retrying (attempt "
                    << retry + 1 << ").");
    outputs = predictVoxels(simulstate, eta, energy);
    if (outputs.empty()) {
      FCS_MSG_ERROR("predictVoxels returned empty outputs on retry; giving up");
      return event_t {};
    }
    retry++;
  }

  // Fill the event structure with the voxel energies
  std::vector<unsigned int> bin_index_vector;
  std::vector<float> E_vector;

  event_t event;

  long unsigned int layer_index = 0;
  long unsigned int layer = m_used_layers.at(layer_index);

  for (long unsigned int voxel_index = 0; voxel_index < m_nVoxels;
       ++voxel_index)
  {
    if (voxel_index == m_layer_boundaries[layer_index + 1]) {
      layer_index = layer_index + 1;
      layer = m_used_layers.at(layer_index);
    }

    const std::string key = std::to_string(voxel_index);
    auto it = outputs.find(key);
    if (it == outputs.end()) {
      FCS_MSG_ERROR("Voxel key '"
                    << key << "' not found in network outputs. "
                    << "Available keys e.g.: " << outputs.begin()->first);
      return event_t {};
    }
    float voxel_energy = it->second;

    if (voxel_energy > 0) {
      if (event.event_data.size() <= layer) {
        event.event_data.resize(layer + 1);
      }
      event.event_data.at(layer).bin_index_vector.push_back(
          voxel_index - m_layer_boundaries[layer_index]);

      // We need energy fractions, not MeV values
      event.event_data.at(layer).E_vector.push_back(voxel_energy / energy);
    }
  }

  return event;
}

VNetworkBase::NetworkOutputs TFCSMLCalorimeterSimulator::predictVoxels() const
{
  // For testing...
  // This function sets the input dimensionality and the number of predicted
  // layers to work with the currently best photon ML simulation model. 382
  // voxels spanned over the presampler, the three EMB layers and the first HCAL
  // layer. This allows for easier testing calls.

  int nEvents = 1;
  int nVoxels = 382;
  int nLayers = 5;

  std::vector<float> eta_vector(nEvents, 2.0);
  std::vector<float> energy_vector(nEvents, 65536.0);
  std::vector<float> z_shape_vector(nEvents * nVoxels, 0.5);
  std::vector<float> z_energy_vector(nEvents * nLayers, 0.5);

  VNetworkBase::NetworkInputs inputs;

  int i = 0;
  for (float eta : eta_vector) {
    inputs["inn_eta_in"].insert(
        std::pair<std::string, double>("variable_" + std::to_string(i), eta));
    i++;
  }

  i = 0;
  for (float energy : energy_vector) {
    inputs["inn_einc_in"].insert(std::pair<std::string, double>(
        "variable_" + std::to_string(i), energy));
    i++;
  }

  i = 0;
  for (float z_shape : z_shape_vector) {
    inputs["cfm_z_shape"].insert(std::pair<std::string, double>(
        "variable_" + std::to_string(i), z_shape));
    i++;
  }

  i = 0;
  for (float z_energy : z_energy_vector) {
    inputs["inn_z_energy"].insert(std::pair<std::string, double>(
        "variable_" + std::to_string(i), z_energy));
    i++;
  }

  FCS_MSG_DEBUG(VNetworkBase::representNetworkInputs(inputs, 1000));

  if (!m_onnx_model) {
    FCS_MSG_ERROR("predictVoxels called before a simulator was loaded");
    return {};
  }

  VNetworkBase::NetworkOutputs outputs = m_onnx_model->compute(inputs);

  FCS_MSG_DEBUG(VNetworkBase::representNetworkOutputs(outputs, 1000));

  return outputs;
}

void TFCSMLCalorimeterSimulator::Print() const
{
  FCS_MSG_INFO("ONNX AICalorimeterSimulator");
}
