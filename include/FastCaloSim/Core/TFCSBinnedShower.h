// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#pragma once

#include <string>
#include <tuple>
#include <vector>

#include <RtypesCore.h>
#include <TMath.h>
#include <hdf5.h>

#include "FastCaloSim/Core/TFCSBinnedShowerBase.h"
#include "FastCaloSim/Core/TFCSLateralShapeParametrizationHitBase.h"
#include "FastCaloSim/Core/TFCSParametrizationBinnedChain.h"
#include "FastCaloSim/Core/TFCSSimulationState.h"
#include "TFCSBinnedShowerBase.h"

class CaloGeo;

class TFCSBinnedShower : public TFCSBinnedShowerBase
{
public:
  struct layer_t
  {
    std::vector<unsigned int> bin_index_vector;
    std::vector<float> E_vector;
  };

  struct event_t
  {
    std::vector<layer_t> event_data;
    float phi_mod {};
    float center_eta {};
    float e_init {};  // Initial energy of the event
  };

  using eventvector_t = std::vector<event_t>;

  struct layer_bins_t
  {
    std::vector<float> R_lower;
    std::vector<float> R_size;
    std::vector<float> alpha_lower;
    std::vector<float> alpha_size;
  };

  using event_bins_t = std::vector<layer_bins_t>;

  TFCSBinnedShower(const char* name = nullptr, const char* title = nullptr);

  ~TFCSBinnedShower() override;

  auto is_match_Ekin_bin(int /*Ekin_bin*/) const -> bool override
  {
    return true;
  };

  void enable_event_matching()
  {
    m_use_event_matching = true;
    m_use_event_cherry_picking = false;
    m_use_eta_matching = false;
  };

  void disable_event_matching() { m_use_event_matching = false; }

  void enable_event_cherry_picking()
  {
    m_use_event_cherry_picking = true;
    m_use_event_matching = false;
    m_use_eta_matching = false;
  };

  void disable_event_cherry_picking() { m_use_event_cherry_picking = false; };

  void enable_eta_matching()
  {
    m_use_eta_matching = true;
    m_use_event_matching = false;
    m_use_event_cherry_picking = false;
  };

  void disable_eta_matching() { m_use_eta_matching = false; };

  void enable_upscaling() { m_use_upscaling = true; };
  void disable_upscaling() { m_use_upscaling = false; };

  // Loads all required data from the given HDF5 file
  // TODO: Define HDF5 file format somewhere
  void load_event_library(const std::string& filename,
                          std::vector<long unsigned int>& layers,
                          bool only_load_meta_data = false);

  void load_sub_bin_distribution(const std::string& filename);

  // If the HDF5 file path is set, the event library will be loaded
  // automatically at the beginning of the simulation.
  // If this class is used, it would be better to use load_event_library()
  // during the param file creation. In this case, the event library
  // is stored and not loaded on the fly.
  void set_hdf5_path(const std::string& filename) { m_hdf5_file = filename; }
  void delete_hdf5_path() { m_hdf5_file.clear(); }
  auto get_hdf5_path() const -> const std::string& { return m_hdf5_file; }

  // Allows to set the layer energy for the given layer and event manually.
  void set_layer_energy(long unsigned int event_index,
                        long unsigned int layer_index,
                        const std::vector<unsigned int>& bin_index_vector,
                        const std::vector<float>& E_vector);

  // Allows to set the voxel boundaries for the given layer manually.
  void set_bin_boundaries(long unsigned int layer_index,
                          const std::vector<float>& R_lower,
                          const std::vector<float>& R_size,
                          const std::vector<float>& alpha_lower,
                          const std::vector<float>& alpha_size);

  // Allows to set the shower center for the given event manually.
  // The layer is needed as reference for the phi modulation calculation.
  // For eta > 1.4 it should be layer 6, otherwise layer 2.
  void set_shower_center_information(long unsigned int event_index,
                                     long unsigned int reference_layer_index,
                                     float eta_center,
                                     float phi_center);

  auto get_eventlibrary() const -> const eventvector_t&
  {
    return m_eventlibrary;
  }

  void set_event_library(const eventvector_t& eventlibrary)
  {
    m_eventlibrary = eventlibrary;
  }

  auto get_coordinates() const -> const event_bins_t& { return m_coordinates; }

  void set_coordinates(const event_bins_t& coordinates)
  {
    m_coordinates = coordinates;
  }

  auto get_sub_bin_distribution() const
      -> const std::vector<std::vector<std::vector<std::vector<float>>>>&
  {
    return m_sub_bin_distribution;
  }

  auto get_upscaling_energies() const -> const std::vector<float>&
  {
    return m_upscaling_energies;
  }

  void set_sub_bin_distribution_and_energies(
      const std::vector<std::vector<std::vector<std::vector<float>>>>&
          sub_bin_distribution,
      const std::vector<float>& upscaling_energies)
  {
    m_sub_bin_distribution = sub_bin_distribution;
    m_upscaling_energies = upscaling_energies;
    m_use_upscaling = true;
  }

protected:
  // Returns the event to be used
  void get_event(TFCSSimulationState& simulstate,
                 float eta_center,
                 float phi_center,
                 float e_init,
                 long unsigned int reference_layer_index) const override;

  // Returns the number of hits that are going to be simulated in the given
  // layer.
  auto get_n_hits(TFCSSimulationState& simulstate,
                  long unsigned int layer_index) const
      -> long unsigned int override;

  // Used to precompute the number of hits for all layers in the event.
  virtual void compute_n_hits_and_elayer(TFCSSimulationState& simulstate) const;

  // Returns the total deposited energy in the given layer for the current event
  auto get_layer_energy(TFCSSimulationState& simulstate,
                        long unsigned int layer_index) const -> float override;

  // Returns the position and energy of the corresponding hit in the given
  // event, layer and bin
  auto get_hit_position_and_energy(TFCSSimulationState& simulstate,
                                   long unsigned int layer_index,
                                   long unsigned int hit_index) const
      -> std::tuple<float, float, float> override;

  // Delete all pointers that were created in get_event()
  void delete_event(TFCSSimulationState& simulstate) const override;

private:
  // Enables to load the event library from an HDF5 file on the fly.
  std::string m_hdf5_file;

  // Store the used event library
  eventvector_t m_eventlibrary;
  event_bins_t m_coordinates;

  // Event matching flag. Ensures that the chose event
  // from the event library is the same as the one taken
  // from the EVNT file. Prevents any bias due to position
  // mismatch. The simulation's event file must be identical
  // to the one used for the event library generation.
  bool m_use_event_matching = false;

  // Event cherry picking flag. Allows select a individual
  // event from the event library based on the shower center.
  // Needed to prevent phi-modulation effects and to artificially
  // shrink the eta slices
  bool m_use_event_cherry_picking = false;

  // Eta matching toggle. Uses the event from the event library that has the
  // most similar eta to the event from the event file. For this setup the event
  // file of the shower extraction and the simulation should NOT be the same!
  bool m_use_eta_matching = false;

  auto find_best_match(float eta_center,
                       float phi_center,
                       float e_init,
                       long unsigned int reference_layer_index,
                       bool phi_mod_matching) const -> long unsigned int;

  auto get_coordinates(TFCSSimulationState& simulstate,
                       long unsigned int layer_index,
                       int bin_index) const -> std::tuple<float, float>;

  // If true, the sub-cell distribution is used to upscale the energy deposition
  bool m_use_upscaling = false;

  // Can be used to interpolate the individual voxels according to a higher
  // resolution average shower
  std::vector<std::vector<std::vector<std::vector<float>>>>
      m_sub_bin_distribution;
  std::vector<float> m_upscaling_energies;  // energies of the avg showers

  void upscale(TFCSSimulationState& simulstate,
               float& R_min,
               float& R_max,
               float& alpha_min,
               float& alpha_max,
               long unsigned int layer_index,
               int bin_index) const;

  auto get_energy_index(TFCSSimulationState& simulstate,
                        long unsigned int layer_index,
                        long unsigned int hit_index) const -> long unsigned int;

  // Helper functions to load the HDF5 dataset
  auto load_hdf5_dataset(const std::string& filename,
                         const std::string& datasetname)
      -> std::tuple<std::vector<float>, std::vector<hsize_t>, bool>;

  void load_layer_energy(const std::string& filename,
                         long unsigned int layer_index);

  void load_bin_boundaries(const std::string& filename,
                           long unsigned int layer_index);

  void load_shower_center_information(const std::string& filename);

  ClassDefOverride(TFCSBinnedShower, 1)  // TFCSBinnedShower
};
