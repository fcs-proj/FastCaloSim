// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#ifndef TFCSBinnedShowerONNX_h
#define TFCSBinnedShowerONNX_h

#include <tuple>
#include <vector>

#include <RtypesCore.h>
#include <TMath.h>
#include <hdf5.h>

#include "FastCaloSim/Core/TFCSBinnedShowerBase.h"
#include "FastCaloSim/Core/TFCSLateralShapeParametrizationHitBase.h"
#include "FastCaloSim/Core/TFCSMLCalorimeterSimulator.h"
#include "FastCaloSim/Core/TFCSParametrizationBinnedChain.h"
#include "FastCaloSim/Core/TFCSSimulationState.h"

class CaloGeo;

class TFCSBinnedShowerONNX : public TFCSBinnedShowerBase
{
public:
  typedef struct
  {
    std::vector<float> R_lower;
    std::vector<float> R_size;
    std::vector<float> alpha_lower;
    std::vector<float> alpha_size;
  } layer_bins_t;

  typedef std::vector<layer_bins_t> event_bins_t;

  TFCSBinnedShowerONNX(const char* name = nullptr, const char* title = nullptr);

  virtual ~TFCSBinnedShowerONNX();

  virtual bool is_match_Ekin_bin(int /*Ekin_bin*/) const override
  {
    return true;
  };

  void enable_upscaling() { m_use_upscaling = true; };
  void disable_upscaling() { m_use_upscaling = false; };

  // Loads the ONNX model for the calorimeter simulation.
  bool load_simulator(const std::string& filename)
  {
    if (m_ai_simulator) {
      delete m_ai_simulator;
      m_ai_simulator = nullptr;
    }
    m_ai_simulator = new TFCSMLCalorimeterSimulator();
    if (!m_ai_simulator->loadSimulator(filename)) {
      delete m_ai_simulator;
      m_ai_simulator = nullptr;
      return false;
    }
    return true;
  }

  // Loads the voxel boundaries and (potentially) the average showers for the
  // upscaling
  // TODO: Define HDF5 file format somewhere
  void load_meta_data(const std::string& filename,
                      std::vector<long unsigned int>& layers);

  void load_sub_bin_distribution(const std::string& filename);

  // Allows to set the voxel boundaries for the given layer manually.
  void set_bin_boundaries(long unsigned int layer_index,
                          std::vector<float> R_lower,
                          std::vector<float> R_size,
                          std::vector<float> alpha_lower,
                          std::vector<float> alpha_size);

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

  // Returns the total deposited energy in the given layer for the current event
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
  // Store the used event library
  event_bins_t m_coordinates;

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

  void load_bin_boundaries(const std::string& filename,
                           long unsigned int layer_index);

  void load_upscaling_distribution(const std::string& filename,
                                   long unsigned int layer_index);

  // Store a reference to the class that handles the actual ONNX call.
  TFCSMLCalorimeterSimulator* m_ai_simulator = nullptr;

  ClassDefOverride(TFCSBinnedShowerONNX, 1)  // TFCSBinnedShowerONNX
};

#endif
