// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#pragma once

// Local includes
#include "FastCaloSim/Core/TFCSLateralShapeParametrizationHitBase.h"

// External includes
#include <string>
#include <tuple>
#include <vector>

#include <RtypesCore.h>

class CaloGeo;

class TFCSPhiModulationCorrection
    : public TFCSLateralShapeParametrizationHitBase
{
public:
  explicit TFCSPhiModulationCorrection(const char* name = nullptr,
                                       const char* title = nullptr);

  ~TFCSPhiModulationCorrection() override;
  TFCSPhiModulationCorrection(const TFCSPhiModulationCorrection& ref) = default;
  TFCSPhiModulationCorrection(TFCSPhiModulationCorrection&& ref) = default;
  auto operator=(const TFCSPhiModulationCorrection& ref)
      -> TFCSPhiModulationCorrection& = default;
  auto operator=(TFCSPhiModulationCorrection&& ref)
      -> TFCSPhiModulationCorrection& = default;

  void load_phi_modulation(const std::string& filename,
                           std::size_t layer_index,
                           float eta_min,
                           float eta_max,
                           float energy_shift = 1.0);

  // Clears only the modulation tables, unlike TFCSParametrization::clear()
  // NOLINTNEXTLINE(bugprone-derived-method-shadowing-base-method): by design
  void clear()
  {
    m_modulation.clear();
    m_min_phi.clear();
    m_min_eta.clear();
    m_energy_shift.clear();
  }

  auto get_phi_modulation() const
      -> const std::vector<std::vector<std::vector<float>>>&
  {
    return m_modulation;
  }
  auto get_min_eta() const -> const std::vector<std::vector<float>>&
  {
    return m_min_eta;
  }
  auto get_energy_shift() const -> const std::vector<std::vector<float>>&
  {
    return m_energy_shift;
  }
  auto get_min_phi() const
      -> const std::vector<std::vector<std::vector<float>>>&
  {
    return m_min_phi;
  }

  void set_phi_modulation_scale(float phi_modulation_scale)
  {
    m_modulation_scale = phi_modulation_scale;
  }

  auto get_phi_modulation_scale() const -> float { return m_modulation_scale; }

  // Adds adds the phi-modulation in the energy to the given hit
  auto add_phi_modulation(Hit& hit) const -> double;

  // Adds adds the phi-modulation in the energy to the given hit
  auto add_phi_modulation(Hit& hit, std::size_t layer_index) const -> double;

  // Adds adds the phi-modulation in the energy to the given hit
  auto add_phi_modulation(double energy,
                          double phi,
                          double eta,
                          std::size_t layer_index) const -> double;

  // Returns the correct eta and phi index for the position of the hit.
  // Used to extract the correct modulation factor from the stored modulation
  // curves.
  auto get_eta_and_phi_index(double phi,
                             double eta,
                             std::size_t layer_index) const
      -> std::tuple<int, std::size_t, std::size_t>;

  void set_geometry(CaloGeo* geo) override
  {
    m_geo = geo;
    TFCSParametrizationBase::set_geometry(geo);
  }

  auto get_geometry() const -> CaloGeo* { return m_geo; }

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
