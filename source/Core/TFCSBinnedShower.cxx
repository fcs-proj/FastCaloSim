// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project
#include <cmath>
#include <numeric>

#include "FastCaloSim/Core/TFCSBinnedShower.h"

#include <H5Cpp.h>

#include "CLHEP/Random/RandFlat.h"
#include "FastCaloSim/Core/TFCSBinnedShowerBase.h"
#include "FastCaloSim/Core/TFCSSimulationState.h"
#include "FastCaloSim/Definitions/FloatCompare.h"
#include "FastCaloSim/Definitions/Phi.h"
#include "FastCaloSim/Geometry/CaloGeo.h"
#include "FastCaloSim/Geometry/LayerRoles.h"
#include "TBuffer.h"
#include "TClass.h"

//=============================================
//======= TFCSBinnedShower =========
//=============================================

TFCSBinnedShower::TFCSBinnedShower(const char* name, const char* title)
    : TFCSBinnedShowerBase(name, title)
{
}

TFCSBinnedShower::~TFCSBinnedShower() = default;

void TFCSBinnedShower::set_layer_energy(
    std::size_t event_index,
    std::size_t layer_index,
    const std::vector<unsigned int>& bin_index_vector,
    const std::vector<float>& E_vector)
{
  // Assert that the event index is valid
  if (event_index >= m_eventlibrary.size()) {
    m_eventlibrary.resize(event_index + 1);
  }

  // Assert that the layer index is valid
  if (layer_index >= m_eventlibrary.at(event_index).event_data.size()) {
    m_eventlibrary.at(event_index).event_data.resize(layer_index + 1);
  }

  // Set the layer energy
  layer_t& layer = m_eventlibrary.at(event_index).event_data.at(layer_index);
  layer.bin_index_vector = bin_index_vector;
  layer.E_vector = E_vector;
}

void TFCSBinnedShower::set_shower_center_information(
    std::size_t event_index,
    std::size_t reference_layer_index,
    double eta_center,
    double phi_center)
{
  // Compute phi_mod
  const double phi_cell_size = m_geo->layer_roles().phi_modulation_period(
      reference_layer_index, eta_center);
  Position pos {};
  pos.m_eta = eta_center;
  pos.m_phi = phi_center;
  const double phi_cell = m_geo->get_cell(reference_layer_index, pos).phi();
  double phi_within_cell = std::fmod(phi_center - phi_cell, phi_cell_size);
  if (phi_within_cell < 0) {
    phi_within_cell += phi_cell_size;
  }

  if (m_eventlibrary.size() <= event_index) {
    m_eventlibrary.resize(event_index + 1);
  }
  if (m_eventlibrary.at(event_index).event_data.size() <= reference_layer_index)
  {
    m_eventlibrary.at(event_index).event_data.resize(reference_layer_index + 1);
  }
  // The event library is stored in single precision
  m_eventlibrary.at(event_index).center_eta = static_cast<float>(eta_center);
  m_eventlibrary.at(event_index).phi_mod = static_cast<float>(phi_within_cell);
}

auto TFCSBinnedShower::find_best_match(double eta_center,
                                       double phi_center,
                                       double e_init,
                                       std::size_t reference_layer_index,
                                       bool phi_mod_matching) const
    -> std::size_t
{
  double phi_cell_size = 0;
  double phi_within_cell = 0;

  if (phi_mod_matching) {
    Position pos {};
    pos.m_eta = eta_center;
    pos.m_phi = phi_center;
    const auto& cell = m_geo->get_cell(reference_layer_index, pos);
    phi_cell_size = m_geo->layer_roles().phi_modulation_period(
        reference_layer_index, eta_center);

    const double phi_cell = cell.phi();
    phi_within_cell = phi_center - phi_cell;

    phi_within_cell = std::fmod(phi_within_cell, phi_cell_size);
    if (phi_within_cell < 0) {
      phi_within_cell += phi_cell_size;
    }
  }

  // Find the event with the closest eta_center and phi_mod (L2 distance)
  double best_distance = std::numeric_limits<double>::max();
  std::size_t best_match = m_eventlibrary.size() + 1;

  for (std::size_t event_index = 0; event_index < m_eventlibrary.size();
       ++event_index)
  {
    // Absolute eta difference. Normalized by the (usual) eta range of 0.05...
    // The library is made for positive eta and used mirrored for negative eta
    const event_t& event = m_eventlibrary.at(event_index);
    const double eta_diff =
        (static_cast<double>(event.center_eta) - std::abs(eta_center)) / 0.05;

    // Logarithmic energy difference. Use power of two (distance between usual
    // energy points) as reasonable distance measure.
    const double e_diff = std::log(static_cast<double>(event.e_init) / e_init)
        / std::numbers::ln2;

    double dist2 = eta_diff * eta_diff + e_diff * e_diff;

    if (phi_mod_matching) {
      // The phases are periodic in phi_cell_size
      const double phi_diff =
          FastCaloSim::periodic_distance(static_cast<double>(event.phi_mod),
                                         phi_within_cell,
                                         phi_cell_size)
          / phi_cell_size;
      dist2 = dist2 + phi_diff * phi_diff;
    }
    const double distance = std::sqrt(dist2);

    if (distance < best_distance) {
      best_distance = distance;
      best_match = event_index;
    }
  }

  if (best_match == m_eventlibrary.size() + 1) {
    FCS_MSG_ERROR("No best match found");
  } else {
    FCS_MSG_INFO("Best match found for eta "
                 << eta_center << " / "
                 << m_eventlibrary.at(best_match).center_eta << " with energy "
                 << e_init << " / " << m_eventlibrary.at(best_match).e_init);
  }

  return best_match;
}

void TFCSBinnedShower::get_event(TFCSSimulationState& simulstate,
                                 double eta_center,
                                 double phi_center,
                                 double e_init,
                                 std::size_t reference_layer_index) const
{
  if (m_eventlibrary.empty()) {
    FCS_MSG_ERROR(
        "No event library loaded. Please load an event library for "
        "TFCSBinnedShower::get_event.");
    return;
  }

  // Uniformly distributed index into the event library
  const auto random_event_index = [&] -> std::size_t
  {
    return static_cast<std::size_t>(std::floor(
        CLHEP::RandFlat::shoot(simulstate.randomEngine(),
                               0,
                               static_cast<double>(m_eventlibrary.size()))));
  };

  const std::size_t event_index = [&] -> std::size_t
  {
    if (simulstate.hasAuxInfo("EventNr"_FCShash) && m_use_event_matching) {
      int const requested = simulstate.getAuxInfo<int>("EventNr"_FCShash);
      if (requested < 0
          || static_cast<std::size_t>(requested) >= m_eventlibrary.size())
      {
        FCS_MSG_WARNING(
            "EventNr " << requested
                       << " is outside the loaded event library of size "
                       << m_eventlibrary.size() << " — falling back to random");
        return random_event_index();
      }
      return static_cast<std::size_t>(requested);
    }
    if (m_use_event_cherry_picking || m_use_eta_matching) {
      // dphi/deta
      const std::size_t best_index =
          find_best_match(eta_center,
                          phi_center,
                          e_init,
                          reference_layer_index,
                          m_use_event_cherry_picking);
      if (best_index >= m_eventlibrary.size()) {
        return random_event_index();
      }
      return best_index;
    }
    return random_event_index();
  }();

  FCS_MSG_DEBUG("Using event index " << event_index << " for eta " << eta_center
                                     << " and phi " << phi_center);

  // Store a copy of the layers of the event
  store_event(simulstate, m_eventlibrary.at(event_index).event_data, e_init);
}

void TFCSBinnedShower::load_event_library(const std::string& filename,
                                          std::vector<std::size_t>& layers,
                                          bool only_load_meta_data)
{
  if (!only_load_meta_data) {
    m_eventlibrary.clear();
  }

  m_coordinates.clear();

  if (m_use_upscaling) {
    m_sub_bin_distribution.clear();
  }

  // layer dependent variables
  for (std::size_t const layer_index : layers) {
    FCS_MSG_INFO("Loading layer " << layer_index << " from file: " << filename);

    // Load the bin boundaries for this layer
    load_bin_boundaries(filename, layer_index);

    if (!only_load_meta_data) {
      // Load the layer energies
      load_layer_energy(filename, layer_index);
    }
  }

  if ((m_use_event_cherry_picking || m_use_eta_matching)
      && !only_load_meta_data)
  {
    load_shower_center_information(filename);
  }
}

void TFCSBinnedShower::load_layer_energy(const std::string& filename,
                                         std::size_t layer_index)
{
  std::string const datasetname = "energy_layer_" + std::to_string(layer_index);

  auto [data, dims, success] = load_hdf5_dataset(filename, datasetname);
  if (!success) {
    FCS_MSG_ERROR("Error while extracting the layer energy for layer "
                  << layer_index << " from " << filename << ".");
    return;
  }

  // Store the data in the event library
  std::vector<unsigned int> bin_index_vector;
  std::vector<float> E_vector;

  for (size_t i = 0; i < data.size(); ++i) {
    std::size_t const event_index = i / dims.at(1);
    const hsize_t bin_index = i % dims.at(1);

    if (bin_index == 0) {
      bin_index_vector.clear();
      E_vector.clear();
    }

    if (!FastCaloSim::exactly_equal(data.at(i), 0.0F)) {
      bin_index_vector.push_back(static_cast<unsigned int>(bin_index));
      E_vector.push_back(data.at(i));
    }

    if (bin_index
        == dims.at(1)
            - 1)  // True for the last voxel of the event in this layer
    {
      set_layer_energy(event_index, layer_index, bin_index_vector, E_vector);
    }
  }
}

void TFCSBinnedShower::load_shower_center_information(
    const std::string& filename)
{
  // Open the HDF5 file
  H5::H5File const file(filename, H5F_ACC_RDONLY);

  // Open the dataset
  std::vector<std::string> datasetnames = {
      "phi_mod", "center_eta", "incident_energy"};

  for (std::size_t i = 0; i < datasetnames.size(); i++) {
    const std::string& datasetname = datasetnames.at(i);
    auto [data, dims, success] = load_hdf5_dataset(filename, datasetname);
    if (!success) {
      if (datasetname == "phi_mod" && m_use_eta_matching) {
        // We do not necessarily need the phi_mod for eta matching, so we can
        // just skip this dataset
        continue;
      }
      FCS_MSG_ERROR("Error while extracting the shower center information from "
                    << filename << "."
                    << "Specifically, the key " << datasetname
                    << " could not be loaded.");

      return;
    }

    for (std::size_t event_index = 0; event_index < data.size(); ++event_index)
    {
      switch (i) {
        case 0:
          m_eventlibrary.at(event_index).phi_mod = data.at(event_index);
          break;
        case 1:
          m_eventlibrary.at(event_index).center_eta = data.at(event_index);
          break;
        case 2:
          m_eventlibrary.at(event_index).e_init = data.at(event_index);
          break;
        default:
          break;
      }
    }
  }
}

void TFCSBinnedShower::set_geometry(CaloGeo* geo)
{
  TFCSBinnedShowerBase::set_geometry(geo);
  if (geo == nullptr || m_hdf5_file.empty() || !m_eventlibrary.empty()) {
    return;
  }
  std::vector<std::size_t> layers(geo->n_layers());
  std::ranges::iota(layers, std::size_t {0});
  FCS_MSG_INFO("Loading event library from " << m_hdf5_file);
  load_event_library(m_hdf5_file, layers);
  FCS_MSG_INFO("Size after loading " << m_eventlibrary.size());
}

void TFCSBinnedShower::Streamer(TBuffer& buffer)
{
  if (buffer.IsReading()) {
    // With an HDF5 file path, the event library is loaded in set_geometry():
    // the geometry is not known yet
    buffer.ReadClassBuffer(TFCSBinnedShower::Class(), this);
    return;
  }
  if (m_hdf5_file.empty()) {
    buffer.WriteClassBuffer(TFCSBinnedShower::Class(), this);
    return;
  }
  // The event library is loaded from the HDF5 file: write the object without
  // it and keep it in memory
  eventvector_t library;
  library.swap(m_eventlibrary);
  buffer.WriteClassBuffer(TFCSBinnedShower::Class(), this);
  library.swap(m_eventlibrary);
}
