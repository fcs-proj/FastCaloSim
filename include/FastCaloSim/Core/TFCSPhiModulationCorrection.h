// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#pragma once

// Local includes
#include "FastCaloSim/Core/TFCSLateralShapeParametrizationHitBase.h"

// External includes
#include <cmath>
#include <string>
#include <tuple>
#include <vector>

#include <RtypesCore.h>
#include <TMath.h>

class CaloGeo;

class TFCSPhiModulationCorrection
    : public TFCSLateralShapeParametrizationHitBase
{
public:
  TFCSPhiModulationCorrection(const char* name = nullptr,
                              const char* title = nullptr);

  ~TFCSPhiModulationCorrection() override;

  void load_phi_modulation(const std::string& filename,
                           long unsigned int layer_index,
                           float eta_min,
                           float eta_max,
                           float energy_shift = 1.0);

  void clear()
  {
    m_modulation.clear();
    m_min_phi.clear();
    m_min_eta.clear();
    m_energy_shift.clear();
  };

  auto get_phi_modulation() const
      -> const std::vector<std::vector<std::vector<float>>>&
  {
    return m_modulation;
  };
  auto get_min_eta() const -> const std::vector<std::vector<float>>&
  {
    return m_min_eta;
  };
  auto get_energy_shift() const -> const std::vector<std::vector<float>>&
  {
    return m_energy_shift;
  };
  auto get_min_phi() const
      -> const std::vector<std::vector<std::vector<float>>>&
  {
    return m_min_phi;
  };

  void set_phi_modulation_scale(float phi_modulation_scale)
  {
    m_modulation_scale = phi_modulation_scale;
  };

  auto get_phi_modulation_scale() const -> float { return m_modulation_scale; };

  static auto get_phi_cell_size(long unsigned int layer, float eta) -> float
  {
    if (layer <= 3) {
      return 2 * TMath::Pi() / 1024;
    } else if (std::abs(eta) < 2.5) {
      return 2 * TMath::Pi() / 768;
    } else {
      return 2 * TMath::Pi() / 256;
    }
  }

  // Adds adds the phi-modulation in the energy to the given hit
  auto add_phi_modulation(Hit& hit) const -> float;

  // Adds adds the phi-modulation in the energy to the given hit
  auto add_phi_modulation(Hit& hit, long unsigned int layer_index) const
      -> float;

  // Adds adds the phi-modulation in the energy to the given hit
  auto add_phi_modulation(float energy,
                          float phi,
                          float eta,
                          long unsigned int layer_index) const -> float;

  // Removes the phi-modulation in the energy from the given hit
  auto remove_phi_modulation(Hit& hit) const -> float;

  // Removes the phi-modulation in the energy from the given hit
  auto remove_phi_modulation(Hit& hit, long unsigned int layer_index) const
      -> float;

  // Removes the phi-modulation in the energy from the given hit
  auto remove_phi_modulation(float energy,
                             float phi,
                             float eta,
                             long unsigned int layer_index) const -> float;

  // Returns the correct eta and phi index for the position of the hit.
  // Used to extract the correct modulation factor from the stored modulation
  // curves.
  auto get_eta_and_phi_index(float phi,
                             float eta,
                             long unsigned int layer_index) const
      -> std::tuple<int, long unsigned int, long unsigned int>;

  void set_geometry(CaloGeo* geo) override
  {
    m_geo = geo;
    TFCSParametrizationBase::set_geometry(geo);
  };

  auto get_geometry() const -> CaloGeo* { return m_geo; };

  auto simulate_hit(Hit& hit,
                    TFCSSimulationState& simulstate,
                    const TFCSTruthState* truth,
                    const TFCSExtrapolationState* extrapol)
      -> FCSReturnCode override;

protected:
  CaloGeo* m_geo = nullptr;  //! do not persistify

private:
  // Needed to reapply the phi modulation
  float m_modulation_scale = 1.0;
  std::vector<std::vector<std::vector<float>>> m_modulation;
  std::vector<std::vector<std::vector<float>>> m_min_phi;
  std::vector<std::vector<float>> m_min_eta;
  std::vector<std::vector<float>> m_energy_shift;

  ClassDefOverride(TFCSPhiModulationCorrection,
                   1)  // TFCSPhiModulationCorrection
};
