// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#include <algorithm>
#include <cmath>
#include <cstdlib>

#include "FastCaloSim/Core/TFCSBinnedShowerBase.h"

#include <CLHEP/Random/RanluxEngine.h>
#include <H5Cpp.h>
#include <TFile.h>
#include <TH2.h>

#include "CLHEP/Random/RandFlat.h"
#include "FastCaloSim/Core/ICaloGeometry.h"
#include "FastCaloSim/Core/TFCSExtrapolationState.h"
#include "FastCaloSim/Core/TFCSLateralShapeParametrizationHitBase.h"
#include "FastCaloSim/Core/TFCSSimulationState.h"
#include "FastCaloSim/Core/TFCSTruthState.h"
#include "HepPDT/ParticleData.hh"
#include "TBuffer.h"
#include "TClass.h"

//=============================================
//======= TFCSBinnedShowerBase =========
//=============================================

TFCSBinnedShowerBase::TFCSBinnedShowerBase(const char* name, const char* title)
    : TFCSLateralShapeParametrizationHitBase(name, title)
{
  reset_OnlyScaleEnergy();
}

TFCSBinnedShowerBase::~TFCSBinnedShowerBase() {}

FCSReturnCode TFCSBinnedShowerBase::simulate(
    TFCSSimulationState& simulstate,
    const TFCSTruthState* truth,
    const TFCSExtrapolationState* extrapol) const
{
  // select a random event from the library
  float eta_center, phi_center;
  long unsigned int reference_layer_index =
      CaloCell_ID_FCS::CaloSample_FCS::EMB2;
  eta_center =
      extrapol->eta(reference_layer_index, TFCSExtrapolationState::SUBPOS_MID);
  if (eta_center > 1.4) {  // Endcap becomes more relevant
    reference_layer_index = CaloCell_ID_FCS::CaloSample_FCS::EME2;
  }
  // TODO: What about the endcap?
  eta_center =
      extrapol->eta(reference_layer_index, TFCSExtrapolationState::SUBPOS_MID);
  phi_center =
      extrapol->phi(reference_layer_index, TFCSExtrapolationState::SUBPOS_MID);

  // Fill the total energy and layer energies into simulstate
  float Einit;
  const float Ekin = truth->Ekin();

  if (OnlyScaleEnergy())
    Einit = simulstate.E();
  else
    Einit = Ekin;

  // Reset the total energy
  simulstate.set_E(0);

  get_event(simulstate, eta_center, phi_center, Einit, reference_layer_index);

  for (long unsigned int layer_index = 0;
       layer_index < CaloCell_ID_FCS::MaxSample;
       ++layer_index)
  {
    float layer_energy = get_layer_energy(simulstate, layer_index);

    // Reset and set the layer energy
    simulstate.set_E(layer_index, 0);
    simulstate.add_E(layer_index, layer_energy);
  }

  if (simulstate.E() > std::numeric_limits<double>::epsilon()) {
    for (int ilayer = 0; ilayer < CaloCell_ID_FCS::MaxSample; ++ilayer) {
      simulstate.set_Efrac(ilayer, simulstate.E(ilayer) / simulstate.E());
    }
  }
  return FCSSuccess;
}

FCSReturnCode TFCSBinnedShowerBase::simulate_hit(
    Hit& hit,
    TFCSSimulationState& simulstate,
    const TFCSTruthState* truth,
    const TFCSExtrapolationState* extrapol)
{
  // Extrapol unused, but needed for the interface
  (void)extrapol;

  const int pdgId = truth->pdgid();
  const float charge = HepPDT::ParticleID(pdgId).charge();
  long unsigned int layer_index = calosample();

  const double center_eta = hit.center_eta();
  const double center_phi = hit.center_phi();
  const double center_r = hit.center_r();
  const double center_z = hit.center_z();

  ATH_MSG_VERBOSE(" Layer " << layer_index << " Extrap eta " << center_eta
                            << " phi " << center_phi << " R " << center_r);

  const float dist000 = TMath::Sqrt(center_r * center_r + center_z * center_z);
  const float eta_jakobi = TMath::Abs(2.0 * TMath::Exp(-center_eta)
                                      / (1.0 + TMath::Exp(-2 * center_eta)));

  long unsigned int hit_index = hit.idx();

  // Get necessary the hit information
  float r, alpha, E;
  std::tie(r, alpha, E) =
      get_hit_position_and_energy(simulstate, layer_index, hit_index);

  hit.reset();
  hit.E() = E;

  if (layer_index <= CaloCell_ID_FCS::CaloSample_FCS::FCAL0) {
    float delta_eta_mm = r * cos(alpha);
    float delta_phi_mm = r * sin(alpha);

    // Particles with negative eta are expected to have the same shape
    // as those with positive eta after transformation: delta_eta -->
    // -delta_eta
    if (center_eta < 0.) {
      delta_eta_mm = -delta_eta_mm;
    }

    // We derive the shower shapes for electrons and positively charged
    // hadrons. Particle with the opposite charge are expected to have the
    // same shower shape after the transformation: delta_phi -->
    // -delta_phi
    if ((charge < 0. && pdgId != 11) || pdgId == -11)
      delta_phi_mm = -delta_phi_mm;

    const float delta_eta = delta_eta_mm / eta_jakobi / dist000;
    const float delta_phi = delta_phi_mm / center_r;

    hit.eta() = center_eta + delta_eta;
    hit.phi() = TVector2::Phi_mpi_pi(center_phi + delta_phi);

    ATH_MSG_VERBOSE(" Hit eta " << hit.eta() << " phi " << hit.phi()
                                << " layer " << layer_index);

  } else {  // FCAL is in (x,y,z)
    const float hit_r = r * cos(alpha) + center_r;
    float delta_phi = r * sin(alpha) / center_r;
    // We derive the shower shapes for electrons and positively charged
    // hadrons. Particle with the opposite charge are expected to have the
    // same shower shape after the transformation: delta_phi -->
    // -delta_phi
    if ((charge < 0. && pdgId != 11) || pdgId == -11)
      delta_phi = -delta_phi;
    const float hit_phi = TVector2::Phi_mpi_pi(center_phi + delta_phi);
    hit.x() = hit_r * cos(hit_phi);
    hit.y() = hit_r * sin(hit_phi);
    hit.z() = center_z;
    ATH_MSG_VERBOSE(" Hit x " << hit.x() << " y " << hit.y() << " layer "
                              << layer_index);
  }

  return FCSSuccess;
}
