// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <memory>
#include <regex>

#include "FastCaloSim/Core/TFCSBinnedShowerBase.h"

#include <H5Cpp.h>
#include <TFile.h>
#include <TKey.h>
#include <TMatrixD.h>

#include "CLHEP/Random/RandFlat.h"
#include "FastCaloSim/Core/TFCSExtrapolationState.h"
#include "FastCaloSim/Core/TFCSLateralShapeParametrizationHitBase.h"
#include "FastCaloSim/Core/TFCSSimulationState.h"
#include "FastCaloSim/Core/TFCSTruthState.h"
#include "FastCaloSim/Definitions/Index.h"
#include "FastCaloSim/Definitions/ParticleData.h"
#include "FastCaloSim/Geometry/CaloGeo.h"

//=============================================
//======= TFCSBinnedShowerBase =========
//=============================================

TFCSBinnedShowerBase::TFCSBinnedShowerBase(const char* name, const char* title)
    : TFCSLateralShapeParametrizationHitBase(name, title)
{
  reset_OnlyScaleEnergy();
}

TFCSBinnedShowerBase::~TFCSBinnedShowerBase() = default;

void TFCSBinnedShowerBase::delete_event(TFCSSimulationState& simulstate)
{
  simulstate.eraseAuxObject("BSEventData"_FCShash);
  simulstate.eraseAuxObject("BSNHits"_FCShash);
  simulstate.eraseAuxObject("BSELayer"_FCShash);
}

auto TFCSBinnedShowerBase::simulate(
    TFCSSimulationState& simulstate,
    const TFCSTruthState* truth,
    const TFCSExtrapolationState* extrapol) const -> FCSReturnCode
{
  // Clean up any event data from a previous simulation call
  delete_event(simulstate);

  // select a random event from the library
  // NOTE: This is ATLAS specific
  // layer 2 is EMB2
  // The event library is made for positive eta; negative eta uses it mirrored
  // (see simulate_hit), so the reference layer is chosen by |eta|
  std::size_t reference_layer = 2;
  if (std::abs(extrapol->eta(reference_layer, FastCaloSim::SubPos::MID)) > 1.4)
  {  // Endcap becomes more relevant
    // NOTE: This is ATLAS specific
    // layer 6 is EME2
    reference_layer = 6;
  }
  // TODO: What about the endcap?
  const double eta_center =
      extrapol->eta(reference_layer, FastCaloSim::SubPos::MID);
  const double phi_center =
      extrapol->phi(reference_layer, FastCaloSim::SubPos::MID);

  // Fill the total energy and layer energies into simulstate
  const double Ekin = truth->Ekin();
  const double Einit = OnlyScaleEnergy() ? simulstate.E() : Ekin;

  // Reset the total energy
  simulstate.set_E(0);

  get_event(simulstate, eta_center, phi_center, Einit, reference_layer);

  for (std::size_t layer_index = 0; layer_index < m_geo->n_layers();
       ++layer_index)
  {
    const double layer_energy = get_layer_energy(simulstate, layer_index);

    // Reset and set the layer energy
    simulstate.set_E(layer_index, 0);
    simulstate.add_E(layer_index, layer_energy);
  }

  if (simulstate.E() > std::numeric_limits<double>::epsilon()) {
    for (std::size_t ilayer = 0; ilayer < m_geo->n_layers(); ++ilayer) {
      simulstate.set_Efrac(ilayer, simulstate.E(ilayer) / simulstate.E());
    }
  }
  return FCSSuccess;
}

auto TFCSBinnedShowerBase::simulate_hit(
    Hit& hit,
    TFCSSimulationState& simulstate,
    const TFCSTruthState* truth,
    const TFCSExtrapolationState* /*extrapol*/) -> FCSReturnCode
{
  const int pdgId = truth->pdgid();
  const double charge = ParticleData::charge(pdgId);
  const std::size_t layer_index = calosample();

  const double center_eta = hit.center_eta();
  const double center_phi = hit.center_phi();
  const double center_r = hit.center_r();
  const double center_z = hit.center_z();

  FCS_MSG_VERBOSE(" Layer " << layer_index << " Extrap eta " << center_eta
                            << " phi " << center_phi << " R " << center_r);

  const double dist000 = FastCaloSim::dist000(center_r, center_z);
  const double eta_jakobi = FastCaloSim::eta_jakobi(center_eta);

  std::size_t const hit_index = hit.idx();

  // Get necessary the hit information
  const auto [r, alpha, E] =
      get_hit_position_and_energy(simulstate, layer_index, hit_index);

  hit.reset();
  hit.set_E(E);

  // NOTE: this is ATLAS dependent
  // layers 0 to 20 are in (eta, phi); the FCAL layers from 21 (FCAL0) on are
  // in (x, y, z), as in TFCSEnergyAndHitGANV2
  if (layer_index <= 20) {
    double delta_eta_mm = r * std::cos(alpha);
    double delta_phi_mm = r * std::sin(alpha);

    // Particles with negative eta are expected to have the same shape
    // as those with positive eta after transformation: delta_eta -->
    // -delta_eta
    if (center_eta < 0.) {
      delta_eta_mm = -delta_eta_mm;
    }

    // We derive the shower shapes for electrons and positively charged
    // hadrons. Particle with the opposite charge are expected to have the
    // same shower shape after the transformation: delta_phi -->
    // -delta_phi
    if ((charge < 0. && pdgId != 11) || pdgId == -11) {
      delta_phi_mm = -delta_phi_mm;
    }

    const double delta_eta = delta_eta_mm / eta_jakobi / dist000;
    const double delta_phi = delta_phi_mm / center_r;

    hit.set_eta_x(center_eta + delta_eta);
    hit.set_phi_y(TVector2::Phi_mpi_pi(center_phi + delta_phi));

    FCS_MSG_VERBOSE(" Hit eta " << hit.eta() << " phi " << hit.phi()
                                << " layer " << layer_index);

  } else {  // FCAL is in (x,y,z)
    const double hit_r = r * std::cos(alpha) + center_r;
    double delta_phi = r * std::sin(alpha) / center_r;
    // We derive the shower shapes for electrons and positively charged
    // hadrons. Particle with the opposite charge are expected to have the
    // same shower shape after the transformation: delta_phi -->
    // -delta_phi
    if ((charge < 0. && pdgId != 11) || pdgId == -11) {
      delta_phi = -delta_phi;
    }
    const double hit_phi = TVector2::Phi_mpi_pi(center_phi + delta_phi);
    hit.set_eta_x(hit_r * std::cos(hit_phi));
    hit.set_phi_y(hit_r * std::sin(hit_phi));
    hit.set_z(center_z);

    FCS_MSG_VERBOSE(" Hit x " << hit.x() << " y " << hit.y() << " layer "
                              << layer_index);
  }

  return FCSSuccess;
}

void TFCSBinnedShowerBase::store_event(TFCSSimulationState& simulstate,
                                       std::vector<layer_t> layers,
                                       double e_init) const
{
  simulstate.setAuxObject("BSEventData"_FCShash, std::move(layers));
  simulstate.setAuxInfo<double>("BSEinit"_FCShash, e_init);
  compute_n_hits_and_elayer(simulstate);
}

void TFCSBinnedShowerBase::set_bin_boundaries(
    std::size_t layer_index,
    const std::vector<float>& R_lower,
    const std::vector<float>& R_size,
    const std::vector<float>& alpha_lower,
    const std::vector<float>& alpha_size)
{
  if (layer_index >= m_coordinates.size()) {
    m_coordinates.resize(layer_index + 1);
  }
  m_coordinates.at(layer_index).R_lower = R_lower;
  m_coordinates.at(layer_index).R_size = R_size;
  m_coordinates.at(layer_index).alpha_lower = alpha_lower;
  m_coordinates.at(layer_index).alpha_size = alpha_size;
}

auto TFCSBinnedShowerBase::get_n_hits(TFCSSimulationState& simulstate,
                                      std::size_t layer_index) const
    -> std::size_t
{
  auto* hits_per_layer_ptr =
      simulstate.getAuxObject<std::vector<std::vector<std::size_t>>>(
          "BSNHits"_FCShash);

  if (hits_per_layer_ptr == nullptr) {
    FCS_MSG_ERROR("Invalid hits per layer information");
    return 0;
  }

  return hits_per_layer_ptr->at(layer_index).back();
}

void TFCSBinnedShowerBase::compute_n_hits_and_elayer(
    TFCSSimulationState& simulstate) const
{
  const auto* layers =
      simulstate.getAuxObject<std::vector<layer_t>>("BSEventData"_FCShash);
  if (layers == nullptr) {
    FCS_MSG_ERROR(
        "compute_n_hits_and_elayer called without BSEventData; "
        "get_event must run first");
    return;
  }
  const auto e_init = simulstate.getAuxInfo<double>("BSEinit"_FCShash);

  // Loop over all layers
  std::size_t const n_layers = layers->size();
  std::vector<std::vector<std::size_t>> hits_per_layer;
  std::vector<double> elayer;
  hits_per_layer.resize(n_layers);
  elayer.resize(n_layers, 0.0);

  for (std::size_t layer_index = 0; layer_index < n_layers; ++layer_index) {
    layer_t const& layer = layers->at(layer_index);

    // Loop over all voxels in the layer
    std::size_t n_hits = 0;
    for (float const e_voxel : layer.E_vector) {
      std::size_t hits_per_bin = 0;
      if (e_voxel > std::numeric_limits<float>::epsilon()) {
        const double e_hits = static_cast<double>(e_voxel) * e_init;
        hits_per_bin = static_cast<std::size_t>(std::min(
            std::max(static_cast<int>(
                         e_hits / static_cast<double>(m_default_hit_energy)),
                     1),
            m_max_hits_per_voxel));
        elayer.at(layer_index) += e_hits;
      }
      n_hits += hits_per_bin;
      hits_per_layer.at(layer_index).push_back(n_hits);
    }
  }
  // Store the hits per layer vector and the energy per layer
  simulstate.setAuxObject("BSNHits"_FCShash, std::move(hits_per_layer));
  simulstate.setAuxObject("BSELayer"_FCShash, std::move(elayer));
}

auto TFCSBinnedShowerBase::get_layer_energy(TFCSSimulationState& simulstate,
                                            std::size_t layer_index) const
    -> double
{
  auto* elayer_ptr =
      simulstate.getAuxObject<std::vector<double>>("BSELayer"_FCShash);
  if (elayer_ptr == nullptr) {
    FCS_MSG_ERROR("Invalid layer energy information");
    return 0.0;
  }
  if (layer_index >= elayer_ptr->size()) {
    return 0.0;
  }

  FCS_MSG_DEBUG("returning energy for layer " << layer_index << ": "
                                              << elayer_ptr->at(layer_index));

  return elayer_ptr->at(layer_index);
}

auto TFCSBinnedShowerBase::get_energy_index(TFCSSimulationState& simulstate,
                                            std::size_t layer_index,
                                            std::size_t hit_index) const
    -> std::size_t
{
  auto* hits_per_layer_ptr =
      simulstate.getAuxObject<std::vector<std::vector<std::size_t>>>(
          "BSNHits"_FCShash);
  if (hits_per_layer_ptr == nullptr) {
    FCS_MSG_ERROR("Invalid hits per layer information");
    return 0;
  }

  if (layer_index >= hits_per_layer_ptr->size()) {
    FCS_MSG_ERROR("Layer index out of bounds: " << layer_index << " >= "
                                                << hits_per_layer_ptr->size());
    return 0;
  }

  // Find the hit index in the hit vector for the given layer
  const std::vector<std::size_t>& hits = hits_per_layer_ptr->at(layer_index);
  auto it = std::ranges::upper_bound(hits, hit_index);
  const auto energy_index = FastCaloSim::index_of(hits, it);

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

auto TFCSBinnedShowerBase::get_hit_position_and_energy(
    TFCSSimulationState& simulstate,
    std::size_t layer_index,
    std::size_t hit_index) const -> std::tuple<double, double, double>
{
  const auto* layers =
      simulstate.getAuxObject<std::vector<layer_t>>("BSEventData"_FCShash);
  if (layers == nullptr) {
    FCS_MSG_ERROR("No BSEventData available in get_hit_position_and_energy");
    return std::make_tuple(0.0, 0.0, 0.0);
  }

  const auto e_init = simulstate.getAuxInfo<double>("BSEinit"_FCShash);

  if (layer_index >= layers->size()) {
    FCS_MSG_ERROR("Layer index out of bounds: " << layer_index
                                                << " >= " << layers->size());
    return std::make_tuple(0.0, 0.0, 0.0);
  }

  std::size_t const energy_index = get_energy_index(
      simulstate, layer_index, hit_index);  // Get the bin index for the hit

  const auto* hits_per_layer_ptr =
      simulstate.getAuxObject<std::vector<std::vector<std::size_t>>>(
          "BSNHits"_FCShash);
  if (hits_per_layer_ptr == nullptr) {
    FCS_MSG_ERROR("No hits per layer available in get_hit_position_and_energy");
    return std::make_tuple(0.0, 0.0, 0.0);
  }

  const std::vector<std::size_t>& hits = hits_per_layer_ptr->at(layer_index);
  const std::size_t hits_per_bin = energy_index == 0
      ? hits.at(energy_index)
      : hits.at(energy_index) - hits.at(energy_index - 1);

  const layer_t& layer = layers->at(layer_index);

  const auto [r, alpha] = get_coordinates(
      simulstate, layer_index, layer.bin_index_vector.at(energy_index));

  const double E = static_cast<double>(layer.E_vector.at(energy_index)) * e_init
      / static_cast<double>(hits_per_bin);

  return std::make_tuple(r, alpha, E);
}

auto TFCSBinnedShowerBase::get_coordinates(TFCSSimulationState& simulstate,
                                           std::size_t layer_index,
                                           unsigned int bin_index) const
    -> std::tuple<double, double>
{
  // The bin coordinates are stored in single precision
  const auto& coordinates = m_coordinates.at(layer_index);
  const auto at_bin = [bin_index](const std::vector<float>& values) -> double
  { return static_cast<double>(values.at(bin_index)); };
  double R_min = at_bin(coordinates.R_lower);
  double R_max = at_bin(coordinates.R_size) + at_bin(coordinates.R_lower);

  double alpha_min = at_bin(coordinates.alpha_lower);
  double alpha_max =
      at_bin(coordinates.alpha_size) + at_bin(coordinates.alpha_lower);

  if (m_use_upscaling) {
    upscale(
        simulstate, R_min, R_max, alpha_min, alpha_max, layer_index, bin_index);
  }
  double R = R_min;  // If the range is too small, just use the minimum value
  if (std::abs(R_max - R_min)
      > static_cast<double>(std::numeric_limits<float>::epsilon()))
  {
    R = CLHEP::RandFlat::shoot(simulstate.randomEngine(), R_min, R_max);
  }
  const double alpha =
      CLHEP::RandFlat::shoot(simulstate.randomEngine(), alpha_min, alpha_max);

  return std::make_tuple(R, alpha);
}

void TFCSBinnedShowerBase::upscale(TFCSSimulationState& simulstate,
                                   double& R_min,
                                   double& R_max,
                                   double& alpha_min,
                                   double& alpha_max,
                                   std::size_t layer_index,
                                   unsigned int bin_index) const
{
  double p = CLHEP::RandFlat::shoot(simulstate.randomEngine(), 0, 1);

  const auto e_init = simulstate.getAuxInfo<double>("BSEinit"_FCShash);
  // Energies and sub-bin probabilities are stored in single precision
  const std::vector<float>& available_energies = m_upscaling_energies;
  const auto energy = [&available_energies](std::size_t index) -> double
  { return static_cast<double>(available_energies.at(index)); };
  const auto sub_bin = [this, layer_index, bin_index](
                           std::size_t index) -> const std::vector<float>&
  { return m_sub_bin_distribution.at(index).at(layer_index).at(bin_index); };

  std::size_t e_index = 0;

  std::array<double, 3> probabilities = {0.25, 0.5, 0.75};
  const auto set_probabilities =
      [&probabilities](const std::vector<float>& sub_bin_probabilities) -> void
  {
    for (std::size_t i = 0; i < probabilities.size(); ++i) {
      probabilities[i] = static_cast<double>(sub_bin_probabilities[i]);
    }
  };

  if (available_energies.size() > 1) {
    // find closest energy index using binary search
    const auto it = std::ranges::upper_bound(
        available_energies,
        e_init,
        std::ranges::less {},
        [](float element) -> double { return static_cast<double>(element); });
    if (it != available_energies.end()) {
      e_index = FastCaloSim::index_of(available_energies, it);
    } else {
      e_index = available_energies.size() - 1;
    }

    const double e_high = energy(e_index);
    if (e_high < e_init || e_index == 0) {
      if (m_sub_bin_distribution.at(e_index).size() > layer_index) {
        set_probabilities(sub_bin(e_index));
      }
    } else {
      if (m_sub_bin_distribution.at(e_index).size() > layer_index
          && m_sub_bin_distribution.at(e_index - 1).size() > layer_index)
      {
        const double e_low = energy(e_index - 1);
        const double f_low =
            std::log(e_high / e_init) / (std::log(e_high / e_low));
        const double f_high = 1 - f_low;
        for (std::size_t i = 0; i < 3; ++i) {
          const auto p_low = static_cast<double>(sub_bin(e_index - 1).at(i));
          const auto p_high = static_cast<double>(sub_bin(e_index).at(i));
          probabilities[i] = f_low * p_low + f_high * p_high;
        }
      }
    }
  }

  else if (available_energies.size() == 1)
  {
    set_probabilities(sub_bin(0));
  }

  const double p_alpha_low =
      probabilities[2] - probabilities[1] + probabilities[0];
  const bool lower_alpha_half = p < p_alpha_low;
  if (lower_alpha_half) {
    alpha_max = (alpha_min + alpha_max) / 2.;
  } else {
    alpha_min = (alpha_min + alpha_max) / 2.;
  }
  double p_r = lower_alpha_half
      ? probabilities[0] / p_alpha_low
      : (probabilities[1] - probabilities[0]) / (1 - p_alpha_low);

  p = CLHEP::RandFlat::shoot(simulstate.randomEngine(), 0, 1);
  if (layer_index != 2) {
    if (p > p_r) {
      R_min = (R_min + R_max) / 2.;
      return;
    }
    R_max = (R_min + R_max) / 2.;
    return;
  }

  // Use linear interpolation for layer 2
  // It works better than uniform sampling for the second layer...
  if (p_r < 0.25) {
    p_r = 0.25;  // Values below 0.25 are not allowed for linear pdf
  } else if (p_r > 0.75) {
    p_r = 0.75;  // Values above 0.75 are not allowed for linear pdf
  } else if (std::abs(p_r - 0.5)
             < static_cast<double>(std::numeric_limits<float>::epsilon()))
  {
    return;  // Best upscaling is uniform sampling. Nothing to do here.
  }

  // Inverse CDF for linear pdf
  double r = (1. - 4. * p_r) / (2. - 4. * p_r)
      - std::sqrt(((1. - 4. * p_r) / (2. - 4. * p_r))
                      * ((1. - 4. * p_r) / (2. - 4. * p_r))
                  + p / ((1. / 2.) - p_r));
  if (r < 0) {
    r = (1. - 4. * p_r) / (2. - 4. * p_r)
        + std::sqrt(((1. - 4. * p_r) / (2. - 4. * p_r))
                        * ((1. - 4. * p_r) / (2. - 4. * p_r))
                    + p / ((1. / 2.) - p_r));
  }

  R_min = R_min + r / 2 * (R_max - R_min);
  R_max = R_min;
}

auto TFCSBinnedShowerBase::load_hdf5_dataset(const std::string& filename,
                                             const std::string& datasetname)
    -> std::tuple<std::vector<float>, std::vector<hsize_t>, bool>
{
  // Open the HDF5 file and dataset
  H5::H5File file(filename, H5F_ACC_RDONLY);

  // check if the dataset exists
  if (!file.exists(datasetname)) {
    return std::make_tuple(
        std::vector<float> {}, std::vector<hsize_t> {}, false);
  }

  H5::DataSet const dataset = file.openDataSet(datasetname);

  // Get the dataspace of the dataset
  H5::DataSpace const dataspace = dataset.getSpace();

  // Get the number of dimensions and the size of each dimension
  int const rank = dataspace.getSimpleExtentNdims();
  std::vector<hsize_t> dims_out(static_cast<std::size_t>(rank));
  dataspace.getSimpleExtentDims(dims_out.data(), nullptr);

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

void TFCSBinnedShowerBase::load_bin_boundaries(const std::string& filename,
                                               std::size_t layer_index)
{
  // Assert that the layer index is valid
  if (layer_index >= m_coordinates.size()) {
    m_coordinates.resize(layer_index + 1);
  }

  std::vector<std::string> datasetnames = {"binstart_radius_layer_",
                                           "binsize_radius_layer_",
                                           "binstart_alpha_layer_",
                                           "binsize_alpha_layer_"};

  for (std::size_t i = 0; i < datasetnames.size(); i++) {
    std::string const datasetname =
        datasetnames.at(i) + std::to_string(layer_index);
    auto [data, dims, success] = load_hdf5_dataset(filename, datasetname);
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
        event_bins.R_lower = std::move(data);
        break;
      case 1:
        event_bins.R_size = std::move(data);
        break;
      case 2:
        event_bins.alpha_lower = std::move(data);
        break;
      case 3:
        event_bins.alpha_size = std::move(data);
        break;
      default:
        break;
    }
  }
}

void TFCSBinnedShowerBase::load_sub_bin_distribution(
    const std::string& filename)
{
  m_use_upscaling = true;
  std::unique_ptr<TFile> file(TFile::Open(filename.c_str(), "READ"));
  if (!file || file->IsZombie()) {
    FCS_MSG_ERROR("Failed to open file: " << filename);
    return;
  }

  std::regex const pattern(R"(probabilities_layer_(\d+)_energy_([0-9.]+))");
  std::map<float, std::map<int, std::vector<std::vector<float>>>> temp_storage;

  TIter next(file->GetListOfKeys());
  while (auto* key = static_cast<TKey*>(next())) {
    std::string const keyname = key->GetName();
    std::smatch match;
    if (std::regex_match(keyname, match, pattern)) {
      int const layer = std::stoi(match[1].str());
      const auto energy = static_cast<float>(std::stod(match[2].str()));

      std::unique_ptr<TMatrixD> matrix(file->Get<TMatrixD>(keyname.c_str()));
      if (matrix) {
        const int n_rows = matrix->GetNrows();
        const int n_cols = matrix->GetNcols();
        std::vector<std::vector<float>> mat_vec;
        mat_vec.reserve(static_cast<std::size_t>(n_rows));
        for (int i = 0; i < n_rows; ++i) {
          std::vector<float> row;
          row.reserve(static_cast<std::size_t>(n_cols));
          for (int j = 0; j < n_cols; ++j) {
            row.push_back(static_cast<float>((*matrix)(i, j)));
          }
          mat_vec.push_back(std::move(row));
        }
        temp_storage[energy][layer] = std::move(mat_vec);
      }
    }
  }

  file->Close();

  // Output containers
  std::vector<float> energies;
  std::vector<std::vector<std::vector<std::vector<float>>>>
      data;  // [energy][layer][row][col]

  for (const auto& [energy, layer_map] : temp_storage) {
    energies.push_back(energy);
    int max_layer = 0;
    for (const auto& [l, _] : layer_map) {
      max_layer = std::max(max_layer, l);
    }

    std::vector<std::vector<std::vector<float>>> layer_vec(
        static_cast<std::size_t>(max_layer) + 1);
    for (const auto& [layer_idx, mat] : layer_map) {
      layer_vec.at(static_cast<std::size_t>(layer_idx)) = mat;
    }
    data.push_back(std::move(layer_vec));
  }

  for (size_t i = 0; i < energies.size(); ++i) {
    FCS_MSG_DEBUG("Energy index " << i << ": " << energies[i] << " GeV");
    for (size_t j = 0; j < data[i].size(); ++j) {
      if (!data[i][j].empty()) {
        FCS_MSG_DEBUG("  Layer " << j << " Shape: (" << data[i][j].size()
                                 << ", " << data[i][j][0].size() << ")");
      }
    }
  }

  m_upscaling_energies = std::move(energies);
  m_sub_bin_distribution = std::move(data);
}
