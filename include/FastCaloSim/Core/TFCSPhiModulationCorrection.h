// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#ifndef TFCSPhiModulationCorrection_h
#define TFCSPhiModulationCorrection_h

// Local includes
#include "FastCaloSim/Core/TFCSLateralShapeParametrizationHitBase.h"

// External includes
#include <fstream>
#include <tuple>
#include <vector>

#include <RtypesCore.h>
#include <TMath.h>

class ICaloGeometry;

class TFCSPhiModulationCorrection
    : public TFCSLateralShapeParametrizationHitBase
{
public:
  TFCSPhiModulationCorrection(const char* name = nullptr,
                              const char* title = nullptr);

  virtual ~TFCSPhiModulationCorrection();

  void load_phi_modulation(std::string filename,
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

  std::vector<std::vector<std::vector<float>>> get_phi_modulation() const
  {
    return m_modulation;
  };
  std::vector<std::vector<float>> get_min_eta() const { return m_min_eta; };
  std::vector<std::vector<float>> get_energy_shift() const
  {
    return m_energy_shift;
  };
  std::vector<std::vector<std::vector<float>>> get_min_phi() const
  {
    return m_min_phi;
  };

  void set_phi_modulation_scale(float phi_modulation_scale)
  {
    m_modulation_scale = phi_modulation_scale;
  };

  float get_phi_modulation_scale() const { return m_modulation_scale; };

  static float get_phi_cell_size(long unsigned int layer, float eta)
  {
    if (layer <= 3) {
      return 2 * TMath::Pi() / 1024;
    } else if (eta < 2.5) {
      return 2 * TMath::Pi() / 768;
    } else {
      return 2 * TMath::Pi() / 256;
    }
  }

  // Adds adds the phi-modulation in the energy to the given hit
  float add_phi_modulation(Hit& hit) const;

  // Adds adds the phi-modulation in the energy to the given hit
  float add_phi_modulation(Hit& hit, long unsigned int layer_index) const;

  // Adds adds the phi-modulation in the energy to the given hit
  float add_phi_modulation(float energy,
                           float phi,
                           float eta,
                           long unsigned int layer_index) const;

  // Removes the phi-modulation in the energy from the given hit
  float remove_phi_modulation(Hit& hit) const;

  // Removes the phi-modulation in the energy from the given hit
  float remove_phi_modulation(Hit& hit, long unsigned int layer_index) const;

  // Removes the phi-modulation in the energy from the given hit
  float remove_phi_modulation(float energy,
                              float phi,
                              float eta,
                              long unsigned int layer_index) const;

  // Returns the correct eta and phi index for the position of the hit.
  // Used to extract the correct modulation factor from the stored modulation
  // curves.
  std::tuple<int, long unsigned int, long unsigned int> get_eta_and_phi_index(
      float phi, float eta, long unsigned int layer_index) const;

  void set_geometry(ICaloGeometry* geo) override
  {
    m_geo = geo;
    TFCSParametrizationBase::set_geometry(geo);
  };

  ICaloGeometry* get_geometry() const { return m_geo; };

  virtual FCSReturnCode simulate_hit(
      Hit& hit,
      TFCSSimulationState& simulstate,
      const TFCSTruthState* truth,
      const TFCSExtrapolationState* extrapol) override;

protected:
  ICaloGeometry* m_geo;  //! do not persistify

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

#endif
