// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#pragma once

#include <string>
#include <tuple>
#include <vector>

#include <RtypesCore.h>
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
  struct event_t
  {
    std::vector<layer_t> event_data;
    float phi_mod {};
    float center_eta {};
    float e_init {};  // Initial energy of the event
  };

  using eventvector_t = std::vector<event_t>;

  explicit TFCSBinnedShower(const char* name = nullptr,
                            const char* title = nullptr);

  ~TFCSBinnedShower() override;
  TFCSBinnedShower(const TFCSBinnedShower&) = default;
  TFCSBinnedShower(TFCSBinnedShower&&) = default;
  auto operator=(const TFCSBinnedShower&) -> TFCSBinnedShower& = default;
  auto operator=(TFCSBinnedShower&&) -> TFCSBinnedShower& = default;

  auto is_match_Ekin_bin(std::size_t /*Ekin_bin*/) const -> bool override
  {
    return true;
  }

  void enable_event_matching()
  {
    m_use_event_matching = true;
    m_use_event_cherry_picking = false;
    m_use_eta_matching = false;
  }

  void disable_event_matching() { m_use_event_matching = false; }

  // Loads the event library from the HDF5 file, if its path is set and the
  // library is not loaded yet, for all layers of the geometry
  void set_geometry(CaloGeo* geo) override;

  void enable_event_cherry_picking()
  {
    m_use_event_cherry_picking = true;
    m_use_event_matching = false;
    m_use_eta_matching = false;
  }

  void disable_event_cherry_picking() { m_use_event_cherry_picking = false; }

  void enable_eta_matching()
  {
    m_use_eta_matching = true;
    m_use_event_matching = false;
    m_use_event_cherry_picking = false;
  }

  void disable_eta_matching() { m_use_eta_matching = false; }

  // Loads all required data from the given HDF5 file
  // TODO: Define HDF5 file format somewhere
  void load_event_library(const std::string& filename,
                          std::vector<std::size_t>& layers,
                          bool only_load_meta_data = false);

  // If the HDF5 file path is set, the event library is not written to a
  // parametrization file; it is loaded from the HDF5 file when the geometry
  // is set (set_geometry(), which is called after reading the file).
  // Alternatively, load_event_library() during the parametrization file
  // creation stores the event library in the file.
  void set_hdf5_path(const std::string& filename) { m_hdf5_file = filename; }
  void delete_hdf5_path() { m_hdf5_file.clear(); }
  auto get_hdf5_path() const -> const std::string& { return m_hdf5_file; }

  // Allows to set the layer energy for the given layer and event manually.
  void set_layer_energy(std::size_t event_index,
                        std::size_t layer_index,
                        const std::vector<unsigned int>& bin_index_vector,
                        const std::vector<float>& E_vector);

  // Allows to set the shower center for the given event manually.
  // The layer is needed as reference for the phi modulation calculation.
  // It should be the layer the simulation takes the shower center from (see
  // FastCaloSim::LayerRoles::shower_center_layers).
  void set_shower_center_information(std::size_t event_index,
                                     std::size_t reference_layer_index,
                                     double eta_center,
                                     double phi_center);

  auto get_eventlibrary() const -> const eventvector_t&
  {
    return m_eventlibrary;
  }

  void set_event_library(const eventvector_t& eventlibrary)
  {
    m_eventlibrary = eventlibrary;
  }

protected:
  // Returns the event to be used
  void get_event(TFCSSimulationState& simulstate,
                 double eta_center,
                 double phi_center,
                 double e_init,
                 std::size_t reference_layer_index) const override;

private:
  // Enables to load the event library from an HDF5 file on the fly.
  std::string m_hdf5_file;

  // Store the used event library
  eventvector_t m_eventlibrary;

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

  auto find_best_match(double eta_center,
                       double phi_center,
                       double e_init,
                       std::size_t reference_layer_index,
                       bool phi_mod_matching) const -> std::size_t;

  void load_layer_energy(const std::string& filename, std::size_t layer_index);

  void load_shower_center_information(const std::string& filename);

  ClassDefOverride(TFCSBinnedShower, 2)  // TFCSBinnedShower
};
