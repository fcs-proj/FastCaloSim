// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#include <algorithm>
#include <cmath>

#include "FastCaloSim/Core/TFCSMLCalorimeterSimulator.h"

#include "CLHEP/Random/RandGauss.h"
#include "FastCaloSim/Core/TFCSNetworkFactory.h"

TFCSMLCalorimeterSimulator::TFCSMLCalorimeterSimulator() = default;

TFCSMLCalorimeterSimulator::~TFCSMLCalorimeterSimulator() = default;

auto TFCSMLCalorimeterSimulator::loadSimulator(const std::string& filename)
    -> bool
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

auto TFCSMLCalorimeterSimulator::predictVoxels(TFCSSimulationState& simulstate,
                                               float eta,
                                               float energy) const
    -> VNetworkBase::NetworkOutputs
{
  // Bring eta into the needed range
  eta = std::abs(eta) * 10;

  // Initialize the energy and eta input vectors
  const auto n_events = static_cast<std::size_t>(m_nEvents);
  std::vector<float> const eta_vector(n_events, eta);
  std::vector<float> const energy_vector(n_events, energy);

  // sample the z vectors according to a standard normal distribution, in the
  // single precision of the network inputs
  std::vector<float> z_shape_vector(n_events * m_nVoxels, 0.0);
  std::vector<float> z_energy_vector(n_events * m_nLayers, 0.0);
  for (auto& z_shape : z_shape_vector) {
    z_shape = static_cast<float>(
        CLHEP::RandGauss::shoot(simulstate.randomEngine(), 0.0, 1.0));
  }
  for (auto& z_energy : z_energy_vector) {
    z_energy = static_cast<float>(
        CLHEP::RandGauss::shoot(simulstate.randomEngine(), 0.0, 1.0));
  }

  // Prepare the inputs for the network
  VNetworkBase::NetworkInputs inputs;

  int i = 0;
  for (float const thisEta : eta_vector) {
    inputs["inn_eta_in"].insert(std::pair<std::string, double>(
        "variable_" + std::to_string(i), thisEta));
    i++;
  }

  i = 0;
  for (float const thisEnergy : energy_vector) {
    inputs["inn_einc_in"].insert(std::pair<std::string, double>(
        "variable_" + std::to_string(i), thisEnergy));
    i++;
  }

  i = 0;
  for (float const z_shape : z_shape_vector) {
    inputs["cfm_z_shape"].insert(std::pair<std::string, double>(
        "variable_" + std::to_string(i), z_shape));
    i++;
  }

  i = 0;
  for (float const z_energy : z_energy_vector) {
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

auto TFCSMLCalorimeterSimulator::getEvent(TFCSSimulationState& simulstate,
                                          float eta,
                                          float energy) const
    -> TFCSMLCalorimeterSimulator::event_t
{
  // Get the voxel energies
  VNetworkBase::NetworkOutputs outputs = predictVoxels(simulstate, eta, energy);

  if (outputs.empty()) {
    FCS_MSG_ERROR(
        "predictVoxels returned empty outputs; returning empty event");
    return event_t {};
  }

  // Check if any output is NaN or infinite; retry up to 5 times if so
  auto has_non_finite = [](const VNetworkBase::NetworkOutputs& out) -> bool
  {
    return std::ranges::any_of(
        out, [](const auto& kv) -> bool { return !std::isfinite(kv.second); });
  };

  int retry = 0;
  while (has_non_finite(outputs)) {
    if (retry >= 5) {
      // NaN voxels would be dropped as voxels without energy and infinite
      // ones kept as infinite energy fractions: return no shower rather than
      // a broken one, as for the other errors
      FCS_MSG_ERROR("Network output still contains NaN or infinity after "
                    << retry << " retries; returning empty event");
      return event_t {};
    }
    FCS_MSG_WARNING(
        "Network output contains NaN or infinity. Retrying (attempt "
        << retry + 1 << ").");
    outputs = predictVoxels(simulstate, eta, energy);
    if (outputs.empty()) {
      FCS_MSG_ERROR("predictVoxels returned empty outputs on retry; giving up");
      return event_t {};
    }
    retry++;
  }

  // Fill the event structure with the voxel energies
  std::vector<unsigned int> const bin_index_vector;
  std::vector<float> const E_vector;

  event_t event;

  std::size_t layer_index = 0;
  std::size_t layer = m_used_layers.at(layer_index);

  for (std::size_t voxel_index = 0; voxel_index < m_nVoxels; ++voxel_index) {
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
    const auto voxel_energy = static_cast<float>(it->second);

    if (voxel_energy > 0) {
      if (event.event_data.size() <= layer) {
        event.event_data.resize(layer + 1);
      }
      event.event_data.at(layer).bin_index_vector.push_back(
          static_cast<unsigned int>(voxel_index
                                    - m_layer_boundaries[layer_index]));

      // We need energy fractions, not MeV values
      event.event_data.at(layer).E_vector.push_back(voxel_energy / energy);
    }
  }

  return event;
}

auto TFCSMLCalorimeterSimulator::predictVoxels() const
    -> VNetworkBase::NetworkOutputs
{
  // For testing...
  // This function sets the input dimensionality and the number of predicted
  // layers to work with the currently best photon ML simulation model. 382
  // voxels spanned over the presampler, the three EMB layers and the first HCAL
  // layer. This allows for easier testing calls.

  const std::size_t nEvents = 1;
  const std::size_t nVoxels = 382;
  const std::size_t nLayers = 5;

  std::vector<float> const eta_vector(nEvents, 2.0);
  std::vector<float> const energy_vector(nEvents, 65536.0);
  std::vector<float> const z_shape_vector(nEvents * nVoxels, 0.5);
  std::vector<float> const z_energy_vector(nEvents * nLayers, 0.5);

  VNetworkBase::NetworkInputs inputs;

  int i = 0;
  for (float const thisEta : eta_vector) {
    inputs["inn_eta_in"].insert(std::pair<std::string, double>(
        "variable_" + std::to_string(i), thisEta));
    i++;
  }

  i = 0;
  for (float const thisEnergy : energy_vector) {
    inputs["inn_einc_in"].insert(std::pair<std::string, double>(
        "variable_" + std::to_string(i), thisEnergy));
    i++;
  }

  i = 0;
  for (float const z_shape : z_shape_vector) {
    inputs["cfm_z_shape"].insert(std::pair<std::string, double>(
        "variable_" + std::to_string(i), z_shape));
    i++;
  }

  i = 0;
  for (float const z_energy : z_energy_vector) {
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
