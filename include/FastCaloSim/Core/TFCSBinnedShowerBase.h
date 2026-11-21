// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#pragma once

// local includes
#include "FastCaloSim/Core/TFCSLateralShapeParametrizationHitBase.h"
#include "FastCaloSim/Core/TFCSSimulationState.h"

// External includes
#include <cassert>
#include <tuple>

#include <RtypesCore.h>

class CaloGeo;

class TFCSBinnedShowerBase : public TFCSLateralShapeParametrizationHitBase
{
public:
  TFCSBinnedShowerBase(const char* name = nullptr, const char* title = nullptr);

  ~TFCSBinnedShowerBase() override;

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
    if (nhits == 0)
      return -1;
    return static_cast<int>(nhits);
  }

protected:
  CaloGeo* m_geo = nullptr;  //! do not persistify

  // Called at the beginning of the simulation to store and or generate the
  // needed shower data for the current event
  virtual void get_event(TFCSSimulationState& simulstate,
                         double eta_center,
                         double phi_center,
                         double e_init,
                         std::size_t reference_layer_index) const = 0;

  // Returns the number of bins that are used in the given layer
  virtual auto get_n_hits(TFCSSimulationState& simulstate,
                          std::size_t layer_index) const -> std::size_t = 0;

  /// Returns the total deposited energy in the given layer for the current
  /// event
  virtual auto get_layer_energy(TFCSSimulationState& simulstate,
                                std::size_t layer_index) const -> double = 0;

  // Returns the position and energy of the corresponding hit in the given
  // event, layer and bin
  virtual auto get_hit_position_and_energy(TFCSSimulationState& simulstate,
                                           std::size_t layer_index,
                                           std::size_t hit_index) const
      -> std::tuple<double, double, double> = 0;

  // Delete all pointers that were created in get_event()
  virtual void delete_event(TFCSSimulationState& simulstate) const = 0;

  // What should be the average energy per hit in the library
  float m_default_hit_energy = 4.;
  // What is the maximum number of hits per voxel (for runtime reasons)
  int m_max_hits_per_voxel = 100;

private:
  ClassDefOverride(TFCSBinnedShowerBase, 1)  // TFCSBinnedShowerBase
};
