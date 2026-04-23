// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#include <cmath>
#include <cstdlib>
#include <regex>

#include "FastCaloSim/Core/TFCSBinnedShowerONNX.h"

#include <CLHEP/Random/RanluxEngine.h>
#include <H5Cpp.h>
#include <TFile.h>
#include <TH2.h>
#include <TKey.h>
#include <TMatrixD.h>

#include "CLHEP/Random/RandFlat.h"
#include "FastCaloSim/Core/TFCSBinnedShowerBase.h"
#include "FastCaloSim/Core/TFCSCenterPositionCalculation.h"
#include "FastCaloSim/Core/TFCSExtrapolationState.h"
#include "FastCaloSim/Core/TFCSInitWithEkin.h"
#include "FastCaloSim/Core/TFCSLateralShapeParametrizationHitBase.h"
#include "FastCaloSim/Core/TFCSMLCalorimeterSimulator.h"
#include "FastCaloSim/Core/TFCSSimulationState.h"
#include "FastCaloSim/Core/TFCSTruthState.h"
#include "FastCaloSim/Definitions/ParticleData.h"
#include "FastCaloSim/Geometry/CaloGeo.h"
#include "TBuffer.h"
#include "TClass.h"

//=============================================
//======= TFCSBinnedShowerONNX =========
//=============================================

// using namespace TFCSBinnedShowerEventTypes;

TFCSBinnedShowerONNX::TFCSBinnedShowerONNX(const char* name, const char* title)
    : TFCSBinnedShowerBase(name, title)
{
}

TFCSBinnedShowerONNX::~TFCSBinnedShowerONNX()
{
  // Delete the AI simulator
  if (m_ai_simulator) {
    delete m_ai_simulator;
    m_ai_simulator = nullptr;
  }
}

void TFCSBinnedShowerONNX::get_event(
    TFCSSimulationState& simulstate,
    float eta_center,
    float phi_center,
    float e_init,
    long unsigned int reference_layer_index) const
{
  (void)reference_layer_index;  // Unused parameter
  (void)phi_center;  // Unused parameter

  const TFCSMLCalorimeterSimulator::event_t& event =
      m_ai_simulator->getEvent(simulstate, eta_center, e_init);

  // Store a pointer to the event
  TFCSMLCalorimeterSimulator::event_t* event_ptr =
      new TFCSMLCalorimeterSimulator::event_t(event);
  simulstate.setAuxInfo<void*>("BSEventData"_FCShash, event_ptr);
  simulstate.setAuxInfo<float>("BSEinit"_FCShash, e_init);

  compute_n_hits_and_elayer(simulstate);
}

void TFCSBinnedShowerONNX::load_meta_data(
    const std::string& filename, std::vector<long unsigned int>& layers)
{
  m_coordinates.clear();

  if (m_use_upscaling) {
    m_sub_bin_distribution.clear();
  }

  // layer dependent variables
  for (long unsigned int layer_index : layers) {
    FCS_MSG_INFO("Loading layer " << layer_index << " from file: " << filename);

    // Load the bin boundaries for this layer
    load_bin_boundaries(filename, layer_index);
  }

  return;
}

std::tuple<std::vector<float>, std::vector<hsize_t>, bool>
TFCSBinnedShowerONNX::load_hdf5_dataset(const std::string& filename,
                                        const std::string& datasetname)
{
  // Open the HDF5 file and dataset
  H5::H5File file(filename, H5F_ACC_RDONLY);

  // check if the dataset exists
  if (!file.exists(datasetname)) {
    return std::make_tuple(
        std::vector<float> {}, std::vector<hsize_t> {}, false);
  }

  H5::DataSet dataset = file.openDataSet(datasetname);

  // Get the dataspace of the dataset
  H5::DataSpace dataspace = dataset.getSpace();

  // Get the number of dimensions and the size of each dimension
  int rank = dataspace.getSimpleExtentNdims();
  std::vector<hsize_t> dims_out(rank);
  dataspace.getSimpleExtentDims(dims_out.data(), NULL);

  // Calculate the total number of elements
  hsize_t totalSize = 1;
  for (const auto& dim : dims_out) {
    totalSize *= dim;
  }

  // Read the dataset into a buffer
  std::vector<float> data(totalSize);
  dataset.read(data.data(), H5::PredType::NATIVE_FLOAT);
  file.close();
  return std::make_tuple(data, dims_out, true);
}

void TFCSBinnedShowerONNX::load_bin_boundaries(const std::string& filename,
                                               long unsigned int layer_index)
{
  // Assert that the layer index is valid
  if (layer_index >= m_coordinates.size()) {
    m_coordinates.resize(layer_index + 1);
  }

  std::vector<std::string> datasetnames = {"binstart_radius_layer_",
                                           "binsize_radius_layer_",
                                           "binstart_alpha_layer_",
                                           "binsize_alpha_layer_"};

  for (long unsigned int i = 0; i < datasetnames.size(); i++) {
    std::string datasetname = datasetnames.at(i) + std::to_string(layer_index);
    std::vector<float> data;
    std::vector<hsize_t> dims;
    bool success;
    std::tie(data, dims, success) = load_hdf5_dataset(filename, datasetname);
    if (!success) {
      FCS_MSG_ERROR("Error while extracting the bin boundaries for layer "
                    << layer_index << " from " << filename << "."
                    << "Specifically, the key " << datasetname
                    << " could not be loaded.");
      return;
    }

    // Fill the corresponding vector in the layer_bins_t structure
    auto& event_bins = m_coordinates.at(layer_index);
    switch (i) {
      case 0:
        event_bins.R_lower = data;
        break;
      case 1:
        event_bins.R_size = data;
        break;
      case 2:
        event_bins.alpha_lower = data;
        break;
      case 3:
        event_bins.alpha_size = data;
        break;
    }
  }
}

// TODO: Could probably be removed in the end
void TFCSBinnedShowerONNX::Streamer(TBuffer& R__b)
{
  // Stream an object of class TFCSBinnedShowerONNX

  if (R__b.IsReading()) {
    R__b.ReadClassBuffer(TFCSBinnedShowerONNX::Class(), this);

  } else {
    R__b.WriteClassBuffer(TFCSBinnedShowerONNX::Class(), this);
  }
}

void TFCSBinnedShowerONNX::set_bin_boundaries(long unsigned int layer_index,
                                              std::vector<float> R_lower,
                                              std::vector<float> R_size,
                                              std::vector<float> alpha_lower,
                                              std::vector<float> alpha_size)
{
  if (layer_index >= m_coordinates.size()) {
    m_coordinates.resize(layer_index + 1);
  }
  m_coordinates.at(layer_index).R_lower = R_lower;
  m_coordinates.at(layer_index).R_size = R_size;
  m_coordinates.at(layer_index).alpha_lower = alpha_lower;
  m_coordinates.at(layer_index).alpha_size = alpha_size;
}

void TFCSBinnedShowerONNX::compute_n_hits_and_elayer(
    TFCSSimulationState& simulstate) const
{
  TFCSMLCalorimeterSimulator::event_t* event =
      static_cast<TFCSMLCalorimeterSimulator::event_t*>(
          simulstate.getAuxInfo<void*>("BSEventData"_FCShash));
  if (!event) {
    FCS_MSG_ERROR(
        "compute_n_hits_and_elayer called without BSEventData; "
        "get_event must run first");
    return;
  }
  float e_init = simulstate.getAuxInfo<float>("BSEinit"_FCShash);

  // Loop over all layers
  long unsigned int n_layers = event->event_data.size();
  std::vector<std::vector<long unsigned int>> hits_per_layer;
  std::vector<float> elayer;
  hits_per_layer.resize(n_layers);
  elayer.resize(n_layers, 0.0f);

  for (long unsigned int layer_index = 0; layer_index < n_layers; ++layer_index)
  {
    TFCSMLCalorimeterSimulator::layer_t& layer =
        event->event_data.at(layer_index);

    // Loop over all voxels in the layer
    long unsigned int n_hits = 0;
    for (float e_voxel : layer.E_vector) {
      long unsigned int hits_per_bin;
      if (e_voxel > std::numeric_limits<float>::epsilon()) {
        hits_per_bin =
            std::min(std::max(int(e_voxel * e_init / m_default_hit_energy), 1),
                     m_max_hits_per_voxel);
        elayer.at(layer_index) += e_voxel * e_init;

      } else {
        hits_per_bin = 0;
      }
      n_hits += hits_per_bin;
      hits_per_layer.at(layer_index).push_back(n_hits);
    }
  }
  // Store the hits per layer vector
  std::vector<std::vector<long unsigned int>>* hits_per_layer_ptr =
      new std::vector<std::vector<long unsigned int>>(hits_per_layer);
  simulstate.setAuxInfo<void*>("BSNHits"_FCShash, hits_per_layer_ptr);

  // Store the energy per layer
  std::vector<float>* elayer_ptr = new std::vector<float>(elayer);
  simulstate.setAuxInfo<void*>("BSELayer"_FCShash, elayer_ptr);
}

long unsigned int TFCSBinnedShowerONNX::get_n_hits(
    TFCSSimulationState& simulstate, long unsigned int layer_index) const
{
  std::vector<std::vector<long unsigned int>>* hits_per_layer_ptr =
      static_cast<std::vector<std::vector<long unsigned int>>*>(
          simulstate.getAuxInfo<void*>("BSNHits"_FCShash));

  if (!hits_per_layer_ptr) {
    FCS_MSG_ERROR("Invalid hits per layer information");
    return 0;
  }

  return hits_per_layer_ptr->at(layer_index).back();
}

float TFCSBinnedShowerONNX::get_layer_energy(
    TFCSSimulationState& simulstate, long unsigned int layer_index) const
{
  std::vector<float>* elayer_ptr = static_cast<std::vector<float>*>(
      simulstate.getAuxInfo<void*>("BSELayer"_FCShash));
  if (!elayer_ptr) {
    FCS_MSG_ERROR("Invalid layer energy information");
    return 0.0f;
  }
  if (layer_index >= elayer_ptr->size()) {
    return 0.0f;
  }

  FCS_MSG_DEBUG("returning energy for layer " << layer_index << ": "
                                              << elayer_ptr->at(layer_index));

  return elayer_ptr->at(layer_index);
}

long unsigned int TFCSBinnedShowerONNX::get_energy_index(
    TFCSSimulationState& simulstate,
    long unsigned int layer_index,
    long unsigned int hit_index) const
{
  std::vector<std::vector<long unsigned int>>* hits_per_layer_ptr =
      static_cast<std::vector<std::vector<long unsigned int>>*>(
          simulstate.getAuxInfo<void*>("BSNHits"_FCShash));
  if (!hits_per_layer_ptr) {
    FCS_MSG_ERROR("Invalid hits per layer information");
    return 0;
  }

  if (layer_index >= hits_per_layer_ptr->size()) {
    FCS_MSG_ERROR("Layer index out of bounds: " << layer_index << " >= "
                                                << hits_per_layer_ptr->size());
    return 0;
  }

  // Find the hit index in the hit vector for the given layer
  const std::vector<long unsigned int>& hits =
      hits_per_layer_ptr->at(layer_index);
  auto it = std::upper_bound(hits.begin(), hits.end(), hit_index);
  long unsigned int energy_index = std::distance(hits.begin(), it);

  if (energy_index >= hits.size()) {
    FCS_MSG_ERROR("Energy index out of bounds: " << energy_index
                                                 << " >= " << hits.size());
    // Print full hits for debugging
    FCS_MSG_ERROR("Hits per layer: ");
    for (const auto& hit : hits) {
      FCS_MSG_ERROR("Hit: " << hit);
    }
    return 0;
  }

  return energy_index;
}

std::tuple<float, float> TFCSBinnedShowerONNX::get_coordinates(
    TFCSSimulationState& simulstate,
    long unsigned int layer_index,
    int bin_index) const
{
  float R_min = m_coordinates.at(layer_index).R_lower.at(bin_index);
  float R_max = m_coordinates.at(layer_index).R_size.at(bin_index)
      + m_coordinates.at(layer_index).R_lower.at(bin_index);

  float alpha_min = m_coordinates.at(layer_index).alpha_lower.at(bin_index);
  float alpha_max = m_coordinates.at(layer_index).alpha_size.at(bin_index)
      + m_coordinates.at(layer_index).alpha_lower.at(bin_index);

  if (m_use_upscaling) {
    upscale(
        simulstate, R_min, R_max, alpha_min, alpha_max, layer_index, bin_index);
  }

  float R = CLHEP::RandFlat::shoot(simulstate.randomEngine(), R_min, R_max);
  float alpha =
      CLHEP::RandFlat::shoot(simulstate.randomEngine(), alpha_min, alpha_max);

  return std::make_tuple(R, alpha);
}

void TFCSBinnedShowerONNX::upscale(TFCSSimulationState& simulstate,
                                   float& R_min,
                                   float& R_max,
                                   float& alpha_min,
                                   float& alpha_max,
                                   long unsigned int layer_index,
                                   int bin_index) const
{
  float p = CLHEP::RandFlat::shoot(simulstate.randomEngine(), 0, 1);

  float e_init = simulstate.getAuxInfo<float>("BSEinit"_FCShash);
  std::vector<float> available_energies = m_upscaling_energies;

  unsigned int e_index = 0;

  std::vector<float> probabilities = {0.25f, 0.5f, 0.75f};

  if (available_energies.size() > 1) {
    // find closest energy index using binary search
    auto it = std::upper_bound(
        available_energies.begin(), available_energies.end(), e_init);
    if (it != available_energies.end()) {
      e_index = std::distance(available_energies.begin(), it);
    } else {
      e_index = available_energies.size() - 1;
    }

    float e_high = available_energies.at(e_index);
    if (e_high < e_init || e_index == 0) {
      if (m_sub_bin_distribution.at(e_index).size() > layer_index) {
        probabilities =
            m_sub_bin_distribution.at(e_index).at(layer_index).at(bin_index);
      }
    } else {
      if (m_sub_bin_distribution.at(e_index).size() > layer_index
          && m_sub_bin_distribution.at(e_index - 1).size() > layer_index)
      {
        float e_low = available_energies.at(e_index - 1);
        float f_low = std::log(e_high / e_init) / (std::log(e_high / e_low));
        float f_high = 1 - f_low;
        for (unsigned int i = 0; i < 3; ++i) {
          float p_low = m_sub_bin_distribution.at(e_index - 1)
                            .at(layer_index)
                            .at(bin_index)
                            .at(i);
          float p_high = m_sub_bin_distribution.at(e_index)
                             .at(layer_index)
                             .at(bin_index)
                             .at(i);
          probabilities[i] = f_low * p_low + f_high * p_high;
        }
      }
    }
  }

  else if (available_energies.size() == 1)
  {
    probabilities = m_sub_bin_distribution.at(0).at(layer_index).at(bin_index);
  }

  float p_alpha_low = probabilities[2] - probabilities[1] + probabilities[0];
  float p_r;

  if (p < p_alpha_low) {
    alpha_max = (alpha_min + alpha_max) / 2.;
    p_r = probabilities[0] / (p_alpha_low);
  } else {
    alpha_min = (alpha_min + alpha_max) / 2.;
    p_r = (probabilities[1] - probabilities[0]) / (1 - p_alpha_low);
  }

  p = CLHEP::RandFlat::shoot(simulstate.randomEngine(), 0, 1);
  if (layer_index != 2) {
    // if ((layer_index != 2) && (layer_index != 0)) {
    if (p > p_r) {
      R_min = (R_min + R_max) / 2.;
      return;
    } else {
      R_max = (R_min + R_max) / 2.;
      return;
    }
  }

  // Use linear interpolation for layer 2
  // It works better than uniform sampling for the second layer...
  if (p_r < 0.25) {
    p_r = 0.25;  // Values below 0.25 are not allowed for linear pdf
  } else if (p_r > 0.75) {
    p_r = 0.75;  // Values above 0.75 are not allowed for linear pdf
  } else if (TMath::Abs(p_r - 0.5) < std::numeric_limits<float>::epsilon()) {
    return;  // Best upscaling is uniform sampling. Nothing to do here.
  }

  // Inverse CDF for linear pdf
  float r = (1. - 4. * p_r) / (2. - 4. * p_r)
      - TMath::Sqrt(((1. - 4. * p_r) / (2. - 4. * p_r))
                        * ((1. - 4. * p_r) / (2. - 4. * p_r))
                    + p / ((1. / 2.) - p_r));
  if (r < 0) {
    r = (1. - 4. * p_r) / (2. - 4. * p_r)
        + TMath::Sqrt(((1. - 4. * p_r) / (2. - 4. * p_r))
                          * ((1. - 4. * p_r) / (2. - 4. * p_r))
                      + p / ((1. / 2.) - p_r));
  }

  R_min = R_min + r / 2 * (R_max - R_min);
  R_max = R_min;
}

std::tuple<float, float, float>
TFCSBinnedShowerONNX::get_hit_position_and_energy(
    TFCSSimulationState& simulstate,
    long unsigned int layer_index,
    long unsigned int hit_index) const
{
  TFCSMLCalorimeterSimulator::event_t* event =
      static_cast<TFCSMLCalorimeterSimulator::event_t*>(
          simulstate.getAuxInfo<void*>("BSEventData"_FCShash));
  if (!event) {
    FCS_MSG_ERROR("No BSEventData available in get_hit_position_and_energy");
    return std::make_tuple(0.0f, 0.0f, 0.0f);
  }

  float e_init = simulstate.getAuxInfo<float>("BSEinit"_FCShash);

  if (layer_index >= event->event_data.size()) {
    FCS_MSG_ERROR("Layer index out of bounds: " << layer_index << " >= "
                                                << event->event_data.size());
    return std::make_tuple(0.0f, 0.0f, 0.0f);
  }

  long unsigned int energy_index = get_energy_index(
      simulstate, layer_index, hit_index);  // Get the bin index for the hit

  std::vector<std::vector<long unsigned int>>* hits_per_layer_ptr =
      static_cast<std::vector<std::vector<long unsigned int>>*>(
          simulstate.getAuxInfo<void*>("BSNHits"_FCShash));

  long unsigned int hits_per_bin;
  if (energy_index == 0) {
    hits_per_bin = hits_per_layer_ptr->at(layer_index).at(energy_index);
  } else {
    hits_per_bin = hits_per_layer_ptr->at(layer_index).at(energy_index)
        - hits_per_layer_ptr->at(layer_index).at(energy_index - 1);
  }

  float r, alpha;

  TFCSMLCalorimeterSimulator::layer_t& layer =
      event->event_data.at(layer_index);

  std::tie(r, alpha) = get_coordinates(
      simulstate, layer_index, layer.bin_index_vector.at(energy_index));

  float E = layer.E_vector.at(energy_index) * e_init / hits_per_bin;

  return std::make_tuple(r, alpha, E);
}

void TFCSBinnedShowerONNX::delete_event(TFCSSimulationState& simulstate) const
{
  // Delete the event data
  void* event_ptr = simulstate.getAuxInfo<void*>("BSEventData"_FCShash);
  if (event_ptr) {
    delete static_cast<TFCSMLCalorimeterSimulator::event_t*>(event_ptr);
    simulstate.setAuxInfo<void*>("BSEventData"_FCShash, nullptr);
  }

  void* n_hits_ptr = simulstate.getAuxInfo<void*>("BSNHits"_FCShash);
  if (n_hits_ptr) {
    delete static_cast<std::vector<std::vector<long unsigned int>>*>(
        n_hits_ptr);
    simulstate.setAuxInfo<void*>("BSNHits"_FCShash, nullptr);
  }

  void* elayer_ptr = simulstate.getAuxInfo<void*>("BSELayer"_FCShash);
  if (elayer_ptr) {
    delete static_cast<std::vector<float>*>(elayer_ptr);
    simulstate.setAuxInfo<void*>("BSELayer"_FCShash, nullptr);
  }
}
void TFCSBinnedShowerONNX::load_sub_bin_distribution(
    const std::string& filename)
{
  m_use_upscaling = true;
  TFile* file = TFile::Open(filename.c_str(), "READ");
  if (!file || file->IsZombie()) {
    std::cerr << "Failed to open file: " << filename << std::endl;
    return;
  }

  std::regex pattern(R"(probabilities_layer_(\d+)_energy_([0-9.]+))");
  std::map<float, std::map<int, std::vector<std::vector<float>>>> temp_storage;

  TIter next(file->GetListOfKeys());
  TKey* key;

  while ((key = (TKey*)next())) {
    std::string keyname = key->GetName();
    std::smatch match;
    if (std::regex_match(keyname, match, pattern)) {
      int layer = std::stoi(match[1].str());
      float energy = std::stod(match[2].str());

      TMatrixD* matrix = dynamic_cast<TMatrixD*>(file->Get(keyname.c_str()));
      if (matrix) {
        std::vector<std::vector<float>> mat_vec(
            matrix->GetNrows(), std::vector<float>(matrix->GetNcols()));
        for (int i = 0; i < matrix->GetNrows(); ++i) {
          for (int j = 0; j < matrix->GetNcols(); ++j) {
            mat_vec[i][j] = static_cast<float>((*matrix)(i, j));
          }
        }
        temp_storage[energy][layer] = std::move(mat_vec);
      }
    }
  }

  file->Close();
  delete file;

  // Output containers
  std::vector<float> energies;
  std::vector<std::vector<std::vector<std::vector<float>>>>
      data;  // [energy][layer][row][col]

  for (const auto& [energy, layer_map] : temp_storage) {
    energies.push_back(energy);
    int max_layer = 0;
    for (const auto& [l, _] : layer_map)
      max_layer = std::max(max_layer, l);

    std::vector<std::vector<std::vector<float>>> layer_vec(max_layer + 1);
    for (const auto& [layer_idx, mat] : layer_map) {
      layer_vec[layer_idx] = mat;
    }
    data.push_back(std::move(layer_vec));
  }

  // Example output
  for (size_t i = 0; i < energies.size(); ++i) {
    std::cout << "Energy index " << i << ": " << energies[i] << " GeV\n";
    for (size_t j = 0; j < data[i].size(); ++j) {
      if (!data[i][j].empty()) {
        std::cout << "  Layer " << j << " Shape: (" << data[i][j].size() << ", "
                  << data[i][j][0].size() << ")\n";
      }
    }
  }

  m_upscaling_energies = energies;
  m_sub_bin_distribution = data;
}
