// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#ifndef TFCSBinnedShower_h
#define TFCSBinnedShower_h

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
  typedef struct
  {
    std::vector<unsigned int> bin_index_vector;
    std::vector<float> E_vector;
  } layer_t;

  typedef struct
  {
    std::vector<layer_t> event_data;
    float phi_mod;
    float center_eta;
    float e_init;  // Initial energy of the event
  } event_t;

  typedef std::vector<event_t> eventvector_t;

  typedef struct
  {
    std::vector<float> R_lower;
    std::vector<float> R_size;
    std::vector<float> alpha_lower;
    std::vector<float> alpha_size;
  } layer_bins_t;

  typedef std::vector<layer_bins_t> event_bins_t;

  TFCSBinnedShower(const char* name = nullptr, const char* title = nullptr);

  virtual ~TFCSBinnedShower();

  virtual bool is_match_Ekin_bin(int /*Ekin_bin*/) const override
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

  void set_default_hit_energy(float energy) { m_default_hit_energy = energy; }
  float get_default_hit_energy() const { return m_default_hit_energy; }
  void set_max_hits_per_voxel(int max_hits) { m_max_hits_per_voxel = max_hits; }
  int get_max_hits_per_voxel() const { return m_max_hits_per_voxel; }

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
  const std::string get_hdf5_path() const { return m_hdf5_file; }

  // Allows to set the layer energy for the given layer and event manually.
  void set_layer_energy(long unsigned int event_index,
                        long unsigned int layer_index,
                        const std::vector<unsigned int> bin_index_vector,
                        const std::vector<float> E_vector);

  // Allows to set the voxel boundaries for the given layer manually.
  void set_bin_boundaries(long unsigned int layer_index,
                          std::vector<float> R_lower,
                          std::vector<float> R_size,
                          std::vector<float> alpha_lower,
                          std::vector<float> alpha_size);

  // Allows to set the shower center for the given event manually.
  // The layer is needed as reference for the phi modulation calculation.
  // For eta > 1.4 it should be layer 6, otherwise layer 2.
  void set_shower_center_information(long unsigned int event_index,
                                     long unsigned int reference_layer_index,
                                     float eta_center,
                                     float phi_center);

  eventvector_t get_eventlibrary() { return m_eventlibrary; }

  void set_event_library(eventvector_t eventlibrary)
  {
    m_eventlibrary = eventlibrary;
  }

  event_bins_t get_coordinates() { return m_coordinates; }

  void set_coordinates(event_bins_t coordinates)
  {
    m_coordinates = coordinates;
  }

  std::vector<std::vector<std::vector<std::vector<float>>>>
  get_sub_bin_distribution() const
  {
    return m_sub_bin_distribution;
  }

  std::vector<float> get_upscaling_energies() const
  {
    return m_upscaling_energies;
  }

  void set_sub_bin_distribution_and_energies(
      std::vector<std::vector<std::vector<std::vector<float>>>>
          sub_bin_distribution,
      std::vector<float> upscaling_energies)
  {
    m_sub_bin_distribution = sub_bin_distribution;
    m_upscaling_energies = upscaling_energies;
    m_use_upscaling = true;
  }

protected:
  // Returns the event to be used
  virtual void get_event(
      TFCSSimulationState& simulstate,
      float eta_center,
      float phi_center,
      float e_init,
      long unsigned int reference_layer_index) const override;

  // Returns the number of hits that are going to be simulated in the given
  // layer.
  virtual long unsigned int get_n_hits(
      TFCSSimulationState& simulstate,
      long unsigned int layer_index) const override;

  // Used to precompute the number of hits for all layers in the event.
  virtual void compute_n_hits_and_elayer(TFCSSimulationState& simulstate) const;

  // Returns the position and energy of the corresponding hit in the given
  // event, layer and bin
  virtual float get_layer_energy(TFCSSimulationState& simulstate,
                                 long unsigned int layer_index) const override;

  // Returns the position and energy of the corresponding hit in the given
  // event, layer and bin
  virtual std::tuple<float, float, float> get_hit_position_and_energy(
      TFCSSimulationState& simulstate,
      long unsigned int layer_index,
      long unsigned int hit_index) const override;

  // Delete all pointers that were created in get_event()
  virtual void delete_event(TFCSSimulationState& simulstate) const override;

private:
  // What should be the average energy per hit in the library
  float m_default_hit_energy = 4.;
  // What is the mamimum number of hits per voxel (for runtime reasons)
  int m_max_hits_per_voxel = 100;

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

  long unsigned int find_best_match(float eta_center,
                                    float phi_center,
                                    float e_init,
                                    long unsigned int reference_layer_index,
                                    bool phi_mod_matching) const;

  std::tuple<float, float> get_coordinates(TFCSSimulationState& simulstate,
                                           long unsigned int layer_index,
                                           int bin_index) const;

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

  long unsigned int get_energy_index(TFCSSimulationState& simulstate,
                                     long unsigned int layer_index,
                                     long unsigned int hit_index) const;

  // Helper functions to load the HDF5 dataset
  std::tuple<std::vector<float>, std::vector<hsize_t>, bool> load_hdf5_dataset(
      const std::string& filename, const std::string& datasetname);

  void load_layer_energy(const std::string& filename,
                         long unsigned int layer_index);

  void load_bin_boundaries(const std::string& filename,
                           long unsigned int layer_index);

  void load_shower_center_information(const std::string& filename);

  ClassDefOverride(TFCSBinnedShower, 1)  // TFCSBinnedShower
};

#endif
