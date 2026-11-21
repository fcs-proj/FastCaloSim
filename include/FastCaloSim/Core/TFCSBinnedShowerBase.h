// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#pragma once

// local includes
#include "FastCaloSim/Core/TFCSLateralShapeParametrizationHitBase.h"
#include "FastCaloSim/Core/TFCSSimulationState.h"

// External includes
#include <cassert>
#include <string>
#include <tuple>
#include <vector>

#include <RtypesCore.h>
#include <hdf5.h>

#include "FastCaloSim/Core/BinnedShowerLayer.h"

class CaloGeo;

class TFCSBinnedShowerBase : public TFCSLateralShapeParametrizationHitBase
{
public:
  /// One layer of a shower: the bins with an energy deposit and their energy
  using layer_t = FastCaloSim::BinnedShowerLayer;

  /// Bin boundaries of one layer, in R and alpha
  struct layer_bins_t
  {
    std::vector<float> R_lower;
    std::vector<float> R_size;
    std::vector<float> alpha_lower;
    std::vector<float> alpha_size;
  };

  using event_bins_t = std::vector<layer_bins_t>;

  explicit TFCSBinnedShowerBase(const char* name = nullptr,
                                const char* title = nullptr);

  ~TFCSBinnedShowerBase() override;
  TFCSBinnedShowerBase(const TFCSBinnedShowerBase&) = default;
  TFCSBinnedShowerBase(TFCSBinnedShowerBase&&) = default;
  auto operator=(const TFCSBinnedShowerBase&)
      -> TFCSBinnedShowerBase& = default;
  auto operator=(TFCSBinnedShowerBase&&) -> TFCSBinnedShowerBase& = default;

  /// Status bit for energy initialization
  enum FCSEnergyInitializationStatusBits
  {
    kOnlyScaleEnergy =
        BIT(18)  ///< Set this bit in the TObject bit field the simulated energy
                 ///< should only be scaled by the GAN
  };

  auto OnlyScaleEnergy() const -> bool { return TestBit(kOnlyScaleEnergy); }

  void set_OnlyScaleEnergy() { SetBit(kOnlyScaleEnergy); }

  void reset_OnlyScaleEnergy() { ResetBit(kOnlyScaleEnergy); }

  void set_default_hit_energy(float energy)
  {
    assert(energy > 0.F && "default hit energy must be > 0");
    m_default_hit_energy = energy;
  }

  void set_max_hits_per_voxel(int max_hits)
  {
    assert(max_hits > 0 && "max hits per voxel must be > 0");
    m_max_hits_per_voxel = max_hits;
  }

  auto get_default_hit_energy() const -> float { return m_default_hit_energy; }

  auto get_max_hits_per_voxel() const -> int { return m_max_hits_per_voxel; }

  // Fill layer energies
  auto simulate(TFCSSimulationState& simulstate,
                const TFCSTruthState* truth,
                const TFCSExtrapolationState* extrapol) const
      -> FCSReturnCode override;

  // Do hit simulation
  auto simulate_hit(Hit& hit,
                    TFCSSimulationState& simulstate,
                    const TFCSTruthState* truth,
                    const TFCSExtrapolationState* extrapol)
      -> FCSReturnCode override;

  void enable_upscaling() { m_use_upscaling = true; }
  void disable_upscaling() { m_use_upscaling = false; }

  void load_sub_bin_distribution(const std::string& filename);

  // Allows to set the voxel boundaries for the given layer manually.
  void set_bin_boundaries(std::size_t layer_index,
                          const std::vector<float>& R_lower,
                          const std::vector<float>& R_size,
                          const std::vector<float>& alpha_lower,
                          const std::vector<float>& alpha_size);

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

  auto get_geometry() -> CaloGeo* { return m_geo; }

  void set_geometry(CaloGeo* geo) override
  {
    m_geo = geo;
    TFCSParametrizationBase::set_geometry(geo);
  }

  auto get_number_of_hits(TFCSSimulationState& simulstate,
                          const TFCSTruthState* /*truth*/,
                          const TFCSExtrapolationState* /*extrapol*/) const
      -> int override
  {
    const std::size_t nhits = get_n_hits(simulstate, calosample());
    if (nhits == 0) {
      return -1;
    }
    return static_cast<int>(nhits);
  }

protected:
  CaloGeo* m_geo = nullptr;  //! do not persistify

  // Called at the beginning of the simulation to select or generate the
  // shower of the current event; stores it with store_event()
  virtual void get_event(TFCSSimulationState& simulstate,
                         double eta_center,
                         double phi_center,
                         double e_init,
                         std::size_t reference_layer_index) const = 0;

  // Stores the layers of the current event and its initial energy in
  // simulstate and computes the number of hits and the energy per layer
  void store_event(TFCSSimulationState& simulstate,
                   std::vector<layer_t> layers,
                   double e_init) const;

  // Returns the number of hits that are simulated in the given layer
  auto get_n_hits(TFCSSimulationState& simulstate,
                  std::size_t layer_index) const -> std::size_t;

  /// Returns the total deposited energy in the given layer for the current
  /// event
  auto get_layer_energy(TFCSSimulationState& simulstate,
                        std::size_t layer_index) const -> double;

  // Returns the position and energy of the corresponding hit in the given
  // event, layer and bin
  auto get_hit_position_and_energy(TFCSSimulationState& simulstate,
                                   std::size_t layer_index,
                                   std::size_t hit_index) const
      -> std::tuple<double, double, double>;

  // Helper functions to load the HDF5 dataset
  static auto load_hdf5_dataset(const std::string& filename,
                                const std::string& datasetname)
      -> std::tuple<std::vector<float>, std::vector<hsize_t>, bool>;

  void load_bin_boundaries(const std::string& filename,
                           std::size_t layer_index);

  // Destroy the data of the previous event that get_event() stored in
  // simulstate
  static void delete_event(TFCSSimulationState& simulstate);

  // What should be the average energy per hit in the library
  float m_default_hit_energy = 4.;
  // What is the maximum number of hits per voxel (for runtime reasons)
  int m_max_hits_per_voxel = 100;

  // Bin boundaries of all layers
  event_bins_t m_coordinates;

  // If true, the sub-cell distribution is used to upscale the energy deposition
  bool m_use_upscaling = false;

  // Can be used to interpolate the individual voxels according to a higher
  // resolution average shower
  std::vector<std::vector<std::vector<std::vector<float>>>>
      m_sub_bin_distribution;
  std::vector<float> m_upscaling_energies;  // energies of the avg showers

private:
  // Precomputes the number of hits and the energy of all layers of the event
  void compute_n_hits_and_elayer(TFCSSimulationState& simulstate) const;

  auto get_coordinates(TFCSSimulationState& simulstate,
                       std::size_t layer_index,
                       unsigned int bin_index) const
      -> std::tuple<double, double>;

  void upscale(TFCSSimulationState& simulstate,
               double& R_min,
               double& R_max,
               double& alpha_min,
               double& alpha_max,
               std::size_t layer_index,
               unsigned int bin_index) const;

  auto get_energy_index(TFCSSimulationState& simulstate,
                        std::size_t layer_index,
                        std::size_t hit_index) const -> std::size_t;

  ClassDefOverride(TFCSBinnedShowerBase, 2)  // TFCSBinnedShowerBase
};
