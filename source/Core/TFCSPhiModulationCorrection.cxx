// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

// Local includes
#include "FastCaloSim/Core/TFCSPhiModulationCorrection.h"

#include "FastCaloSim/Core/TFCSExtrapolationState.h"
#include "FastCaloSim/Core/TFCSSimulationState.h"
#include "FastCaloSim/Core/TFCSTruthState.h"
#include "FastCaloSim/Geometry/CaloGeo.h"

// External includes
#include <RtypesCore.h>
#include <TFile.h>
#include <TH2F.h>
#include <TMath.h>
#include <TParameter.h>

//=============================================
//======= TFCSPhiModulationCorrection =========
//=============================================

TFCSPhiModulationCorrection::TFCSPhiModulationCorrection(const char* name,
                                                         const char* title)
    : TFCSLateralShapeParametrizationHitBase(name, title)
{
}

TFCSPhiModulationCorrection::~TFCSPhiModulationCorrection() {}

void TFCSPhiModulationCorrection::load_phi_modulation(
    std::string filename,
    long unsigned int layer_index,
    float eta_min,
    float eta_max,
    float energy_shift)
{
  if (m_min_eta.size() <= layer_index) {
    m_min_eta.resize(layer_index + 1);
  }
  if (m_min_phi.size() <= layer_index) {
    m_min_phi.resize(layer_index + 1);
  }
  if (m_modulation.size() <= layer_index) {
    m_modulation.resize(layer_index + 1);
  }
  if (m_energy_shift.size() <= layer_index) {
    m_energy_shift.resize(layer_index + 1);
  }
  std::vector<float>& eta_mins = m_min_eta.at(layer_index);
  std::vector<std::vector<float>>& phi_mins = m_min_phi.at(layer_index);
  std::vector<std::vector<float>>& modulation = m_modulation.at(layer_index);
  std::vector<float>& energy_shifts = m_energy_shift.at(layer_index);

  if (eta_mins.empty()) {
    eta_mins.push_back(eta_min);
    eta_mins.push_back(eta_max);
    phi_mins.resize(2);
    modulation.resize(2);
    energy_shifts.push_back(-1.0);
    energy_shifts.push_back(-1.0);
  }

  if (eta_mins.at(0) > eta_min) {
    eta_mins.insert(eta_mins.begin(), eta_min);
    phi_mins.insert(phi_mins.begin(), std::vector<float>());
    modulation.insert(modulation.begin(), std::vector<float>());
    energy_shifts.insert(energy_shifts.begin(), -1.0);
  }

  // Find the matching eta index
  auto eta_it = std::upper_bound(m_min_eta.at(layer_index).begin(),
                                 m_min_eta.at(layer_index).end(),
                                 eta_min);
  long unsigned int eta_index =
      std::distance(m_min_eta.at(layer_index).begin(), eta_it) - 1;

  // Check if we have to append further boundaries
  if (eta_index >= eta_mins.size() - 1) {
    if (eta_min > eta_mins.at(eta_mins.size() - 1)) {
      eta_mins.push_back(eta_min);
      phi_mins.push_back(std::vector<float>());
      modulation.push_back(std::vector<float>());
      energy_shifts.push_back(-1.0);
      eta_index++;
    }
    eta_mins.push_back(eta_max);
    phi_mins.push_back(std::vector<float>());
    modulation.push_back(std::vector<float>());
    energy_shifts.push_back(-1.0);
  } else {
    // Check if we need to insert a new eta_min
    if (eta_mins.at(eta_index) < eta_min) {
      eta_mins.insert(eta_mins.begin() + eta_index + 1, eta_min);
      phi_mins.insert(phi_mins.begin() + eta_index + 1, std::vector<float>());
      modulation.insert(modulation.begin() + eta_index + 1,
                        std::vector<float>());
      energy_shifts.insert(energy_shifts.begin() + eta_index + 1, -1.0);
      eta_index++;
    }
    // Check if we need to insert a new eta_max
    if (eta_mins.at(eta_index + 1) > eta_max) {
      eta_mins.insert(eta_mins.begin() + eta_index + 1, eta_max);
      phi_mins.insert(phi_mins.begin() + eta_index + 1, std::vector<float>());
      modulation.insert(modulation.begin() + eta_index + 1,
                        std::vector<float>());
      energy_shifts.insert(energy_shifts.begin() + eta_index + 1, -1.0);
    }
  }

  energy_shifts.at(eta_index) = energy_shift;

  FCS_MSG_DEBUG("Loading phi modulation correction from " << filename);

  std::unique_ptr<TFile> muon_corr(TFile::Open(filename.c_str()));
  if (!muon_corr || muon_corr->IsZombie()) {
    FCS_MSG_ERROR("Failed to open phi modulation file: " << filename);
    return;
  }

  TH2F* muon_corr_hist = dynamic_cast<TH2F*>(muon_corr->Get("hWt_Layer0"));
  if (!muon_corr_hist) {
    FCS_MSG_ERROR("Failed to load hWt_Layer0 histogram from " << filename);
    return;
  }

  int n_bins = muon_corr_hist->GetNbinsX();

  for (int i = 1; i <= n_bins; i++) {
    phi_mins.at(eta_index).push_back(
        muon_corr_hist->GetXaxis()->GetBinLowEdge(i));
  }
  phi_mins.at(eta_index).push_back(
      muon_corr_hist->GetXaxis()->GetBinUpEdge(n_bins));

  for (int i = 1; i <= n_bins; i++) {
    modulation.at(eta_index).push_back(muon_corr_hist->GetBinContent(i));
  }

  TParameter<double>* param =
      dynamic_cast<TParameter<double>*>(muon_corr->Get("energy_shift"));
  if (param) {
    energy_shifts.at(eta_index) = param->GetVal();
    FCS_MSG_DEBUG("Energy shift: " << energy_shifts.at(eta_index));
  }

  return;
}

// TODO: Use FCSReturnCode
std::tuple<int, long unsigned int, long unsigned int>
TFCSPhiModulationCorrection::get_eta_and_phi_index(
    float phi, float eta, long unsigned int layer_index) const
{
  float eta_abs = TMath::Abs(eta);

  auto eta_it = std::upper_bound(m_min_eta.at(layer_index).begin(),
                                 m_min_eta.at(layer_index).end(),
                                 eta_abs);
  long unsigned int eta_index =
      std::distance(m_min_eta.at(layer_index).begin(), eta_it) - 1;

  if (eta_index >= m_min_eta.at(layer_index).size() - 1) {
    FCS_MSG_ERROR("Found invalid eta index for phi modulation");
    FCS_MSG_ERROR("Layer: " << layer_index);
    FCS_MSG_ERROR("Eta: " << eta);
    FCS_MSG_ERROR("Min eta:" << m_min_eta.at(layer_index).at(0) << " Max eta: "
                             << m_min_eta.at(layer_index).back());
    FCS_MSG_ERROR("Eta index: " << eta_index);
    FCS_MSG_ERROR("Number of eta bins: " << m_min_eta.at(layer_index).size());
    FCS_MSG_ERROR("Eta bin boundaries: ");
    for (const auto& eta_min : m_min_eta.at(layer_index)) {
      FCS_MSG_ERROR("  " << eta_min);
    }
    return std::make_tuple(1, 0, 0);  // Error code 1: Invalid eta index
  }

  if (m_min_eta.at(layer_index).at(eta_index) > eta_abs
      || m_min_eta.at(layer_index).at(eta_index + 1) < eta_abs)
  {
    FCS_MSG_ERROR(
        "Found eta outside of the specified eta range for the phi modulation");
    FCS_MSG_ERROR("Layer: " << layer_index);
    FCS_MSG_ERROR("Eta: " << eta);
    FCS_MSG_ERROR("Phi: " << phi);
    FCS_MSG_ERROR("Eta min of bin: "
                  << m_min_eta.at(layer_index).at(eta_index)
                  << " Eta max of bin: "
                  << m_min_eta.at(layer_index).at(eta_index + 1));
    return std::make_tuple(2, 0, 0);  // Error code 2: Eta outside of range
  }

  Position pos {};
  pos.m_eta = eta;
  pos.m_phi = phi;
  const auto& cell = m_geo->get_cell(layer_index, pos);

  float cell_phi = cell.phi();

  float phi_within_cell = phi - cell_phi;

  float phi_cell_size = get_phi_cell_size(layer_index, eta);

  phi_within_cell = fmod(phi_within_cell, phi_cell_size);

  if (phi_within_cell < 0)
    phi_within_cell += phi_cell_size;

  const std::vector<float>& phi_mins = m_min_phi.at(layer_index).at(eta_index);
  const std::vector<float>& modulation =
      m_modulation.at(layer_index).at(eta_index);

  if (phi_mins.size() == 0) {
    return std::make_tuple(-1, 0, 0);  // Empty phi boundaries -> nothing to do
  }

  // Find the phi bin of the hit
  auto phi_it =
      std::upper_bound(phi_mins.begin(), phi_mins.end(), phi_within_cell);
  long unsigned int phi_index = std::distance(phi_mins.begin(), phi_it) - 1;

  if (phi_index >= modulation.size()) {
    FCS_MSG_ERROR("Found bin "
                  << phi_index
                  << " outside of the modulation correction vector");
    FCS_MSG_ERROR("Phi: " << phi << " Cell phi: " << cell_phi
                          << " Phi within cell: " << phi_within_cell);
    FCS_MSG_ERROR("Eta: " << eta);
    FCS_MSG_ERROR("Phi cell size: " << phi_cell_size);
    FCS_MSG_ERROR("Modulation correction size: " << modulation.size());
    FCS_MSG_ERROR("Phi bin boundaries size: " << phi_mins.size());
    FCS_MSG_ERROR("Last Phi bin boundary: " << phi_mins.back());

    return std::make_tuple(
        3, 0, 0);  // Error code 3: Phi index outside of modulation vector
  }

  return std::make_tuple(0, eta_index, phi_index);  // No error
}

float TFCSPhiModulationCorrection::add_phi_modulation(
    TFCSLateralShapeParametrizationHitBase::Hit& hit) const
{
  return add_phi_modulation(hit, calosample());
}

float TFCSPhiModulationCorrection::add_phi_modulation(
    TFCSLateralShapeParametrizationHitBase::Hit& hit,
    long unsigned int layer_index) const
{
  if (layer_index >= m_min_eta.size()) {
    return hit.E();
  }

  if (m_min_eta.at(layer_index).size() == 0) {
    return hit.E();
  }

  if (m_modulation_scale < std::numeric_limits<float>::epsilon()) {
    return hit.E();
  }

  float energy = hit.E();
  float phi = hit.phi();
  float eta = hit.eta();

  float reweighted_energy = add_phi_modulation(energy, phi, eta, layer_index);
  hit.set_E(reweighted_energy);

  return reweighted_energy;
}

float TFCSPhiModulationCorrection::add_phi_modulation(
    float energy, float phi, float eta, long unsigned int layer_index) const
{
  if (layer_index >= m_min_eta.size()) {
    return energy;
  }

  if (m_min_eta.at(layer_index).size() == 0) {
    return energy;
  }

  if (m_modulation_scale < std::numeric_limits<float>::epsilon()) {
    return energy;
  }

  long unsigned int eta_index;
  long unsigned int phi_index;
  int error_code;
  std::tie(error_code, eta_index, phi_index) =
      get_eta_and_phi_index(phi, eta, layer_index);

  if (error_code != 0) {
    return energy;  // Return original energy in case of error --- no modulation
                    // applied
  }

  const std::vector<float>& modulation =
      m_modulation.at(layer_index).at(eta_index);

  energy *= ((modulation.at(phi_index) - 1) * m_modulation_scale + 1);
  energy *= m_energy_shift.at(layer_index).at(eta_index);

  return energy;
}

float TFCSPhiModulationCorrection::remove_phi_modulation(
    TFCSLateralShapeParametrizationHitBase::Hit& hit) const
{
  return remove_phi_modulation(hit, calosample());
}

float TFCSPhiModulationCorrection::remove_phi_modulation(
    TFCSLateralShapeParametrizationHitBase::Hit& hit,
    long unsigned int layer_index) const
{
  if (layer_index >= m_min_eta.size()) {
    return hit.E();
  }

  if (m_min_eta.at(layer_index).size() == 0) {
    return hit.E();
  }

  if (m_modulation_scale < std::numeric_limits<float>::epsilon()) {
    return hit.E();
  }

  float energy = hit.E();
  float phi = hit.phi();
  float eta = hit.eta();

  float reweighted_energy =
      remove_phi_modulation(energy, phi, eta, layer_index);
  hit.set_E(reweighted_energy);

  return reweighted_energy;
}

float TFCSPhiModulationCorrection::remove_phi_modulation(
    float energy, float phi, float eta, long unsigned int layer_index) const
{
  if (layer_index >= m_min_eta.size()) {
    return energy;
  }

  if (m_min_eta.at(layer_index).size() == 0) {
    return energy;
  }

  if (m_modulation_scale < std::numeric_limits<float>::epsilon()) {
    return energy;
  }

  long unsigned int eta_index;
  long unsigned int phi_index;
  int error_code;
  std::tie(error_code, eta_index, phi_index) =
      get_eta_and_phi_index(phi, eta, layer_index);

  if (error_code != 0) {
    return energy;  // Return original energy in case of error --- no modulation
                    // removed
  }

  const std::vector<float>& modulation =
      m_modulation.at(layer_index).at(eta_index);

  energy /= ((modulation.at(phi_index) - 1) * m_modulation_scale + 1);

  return energy;
}

FCSReturnCode TFCSPhiModulationCorrection::simulate_hit(
    Hit& hit,
    TFCSSimulationState& simulstate,
    const TFCSTruthState* truth,
    const TFCSExtrapolationState* extrapol)
{
  // Extrapol unused, but needed for the interface
  (void)extrapol;

  // Simulstate unused, but needed for the interface
  (void)simulstate;

  // truth unused, but needed for the interface
  (void)truth;

  add_phi_modulation(hit);

  return FCSSuccess;
}
