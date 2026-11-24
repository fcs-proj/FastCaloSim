// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <ranges>
#include <utility>
#include <vector>

#include "FastCaloSim/Core/TFCSHistoLateralShapeParametrization.h"

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandPoisson.h"
#include "FastCaloSim/Core/TFCSExtrapolationState.h"
#include "FastCaloSim/Core/TFCSSimulationState.h"
#include "FastCaloSim/Definitions/Index.h"
#include "FastCaloSim/Definitions/ParticleData.h"
#include "TFile.h"
#include "TH2.h"

//=============================================
//======= TFCSHistoLateralShapeParametrization =========
//=============================================

TFCSHistoLateralShapeParametrization::TFCSHistoLateralShapeParametrization(
    const char* name, const char* title)
    : TFCSLateralShapeParametrizationHitBase(name, title)

{
  TFCSHistoLateralShapeParametrization::reset_phi_symmetric();
}

TFCSHistoLateralShapeParametrization::~TFCSHistoLateralShapeParametrization() =
    default;

void TFCSHistoLateralShapeParametrization::set_geometry(CaloGeo* geo)
{
  TFCSLateralShapeParametrizationHitBase::set_geometry(geo);
  const auto is_nan = [](float c) -> bool { return std::isnan(c); };
  if (std::ranges::none_of(m_hist.get_HistoContents(), is_nan)) {
    return;
  }

  // The contents are a cumulative distribution, which the sampling searches
  // with upper_bound: the repair keeps it non-decreasing. Fixed on a copy,
  // which then replaces the contents of the histogram (and rebuilds its search
  // table).
  std::vector<float> contents = m_hist.get_HistoContents();
  const auto first_valid = std::ranges::find_if_not(contents, is_nan);
  if (first_valid == contents.end()) {
    std::ranges::fill(contents, 1.0F);
    FCS_MSG_WARNING("nan in histo content for "
                    << GetTitle()
                    << " for all bins. Fixed to probability 1 causing hits "
                       "to be deposited in the shower center");
  } else {
    const auto last_valid =
        std::ranges::find_if_not(contents | std::views::reverse, is_nan).base();
    // Before the first valid bin no probability, after the last valid bin all
    // of it, and a bin in between adds none to the bin before it
    std::fill(contents.begin(), first_valid, 0.0F);
    std::fill(last_valid, contents.end(), 1.0F);
    std::size_t n_inner = 0;
    for (auto it = first_valid; it != last_valid; ++it) {
      if (is_nan(*it)) {
        *it = *(it - 1);
        ++n_inner;
      }
    }
    FCS_MSG_WARNING("nan in histo content for "
                    << GetTitle() << ". Fixed "
                    << FastCaloSim::index_of(contents, first_valid)
                    << " bins at the start with probability 0, "
                    << static_cast<std::size_t>(contents.end() - last_valid)
                    << " bins at the end with probability 1 and " << n_inner
                    << " bins in between with the probability of the bin "
                       "before.");
  }
  m_hist.set_HistoContents(std::move(contents));
}

auto TFCSHistoLateralShapeParametrization::get_sigma2_fluctuation(
    TFCSSimulationState& /*simulstate*/,
    const TFCSTruthState* /*truth*/,
    const TFCSExtrapolationState* /*extrapol*/) const -> double
{
  // Limit to factor 1000 fluctuations
  const auto nhits = static_cast<double>(m_nhits);
  if (nhits < 0.001) {
    return 1000;
  }
  return 1.0 / nhits;
}

auto TFCSHistoLateralShapeParametrization::get_number_of_hits(
    TFCSSimulationState& simulstate,
    const TFCSTruthState* /*truth*/,
    const TFCSExtrapolationState* /*extrapol*/) const -> int
{
  if (simulstate.randomEngine() == nullptr) {
    return -1;
  }

  return static_cast<int>(CLHEP::RandPoisson::shoot(
      simulstate.randomEngine(), static_cast<double>(m_nhits)));
}

void TFCSHistoLateralShapeParametrization::set_number_of_hits(float nhits)
{
  m_nhits = nhits;
}

auto TFCSHistoLateralShapeParametrization::simulate_hit(
    Hit& hit,
    TFCSSimulationState& simulstate,
    const TFCSTruthState* truth,
    const TFCSExtrapolationState* /*extrapol*/) -> FCSReturnCode
{
  if (simulstate.randomEngine() == nullptr) {
    return FCSFatal;
  }

  const int pdgId = truth->pdgid();
  const double charge = ParticleData::charge(pdgId);

  const std::size_t cs = calosample();
  const double center_eta = hit.center_eta();
  const double center_phi = hit.center_phi();
  const double center_r = hit.center_r();
  const double center_z = hit.center_z();

  if (std::isnan(center_r) or std::isnan(center_z) or std::isnan(center_eta)
      or std::isnan(center_phi))
  {  // Check if extrapolation fails
    return FCSFatal;
  }

  // Set by m_hist.rnd_to_fct()
  double alpha = 0;
  double r = 0;
  const double rnd0 = CLHEP::RandFlat::shoot(simulstate.randomEngine());
  double rnd1 = CLHEP::RandFlat::shoot(simulstate.randomEngine());
  if (is_phi_symmetric()) {
    if (rnd1 >= 0.5) {  // Fill negative phi half of shape
      rnd1 -= 0.5;
      rnd1 *= 2;
      m_hist.rnd_to_fct(alpha, r, rnd0, rnd1);
      alpha = -alpha;
    } else {  // Fill positive phi half of shape
      rnd1 *= 2;
      m_hist.rnd_to_fct(alpha, r, rnd0, rnd1);
    }
  } else {
    m_hist.rnd_to_fct(alpha, r, rnd0, rnd1);
  }
  if (std::isnan(alpha) || std::isnan(r)) {
    FCS_MSG_ERROR("  Histogram: "
                  << m_hist.get_HistoBordersx().size() - 1 << "*"
                  << m_hist.get_HistoBordersy().size() - 1
                  << " bins, #hits=" << m_nhits << " alpha=" << alpha
                  << " r=" << r << " rnd0=" << rnd0 << " rnd1=" << rnd1);

    FCS_MSG_ERROR("  This error could probably be retried");
    return FCSFatal;
  }

  r *= static_cast<double>(m_r_scale);
  r += static_cast<double>(m_r_offset);
  r = std::max<double>(r, 0);

  double delta_eta_mm = r * std::cos(alpha);
  double delta_phi_mm = r * std::sin(alpha);

  // Particles with negative eta are expected to have the same shape as those
  // with positive eta after transformation: delta_eta --> -delta_eta
  if (center_eta < 0.) {
    delta_eta_mm = -delta_eta_mm;
  }
  // We derive the shower shapes for electrons and positively charged hadrons.
  // Particle with the opposite charge are expected to have the same shower
  // shape after the transformation: delta_phi --> -delta_phi
  if ((charge < 0. && pdgId != 11) || pdgId == -11) {
    delta_phi_mm = -delta_phi_mm;
  }

  const double dist000 = FastCaloSim::dist000(center_r, center_z);
  const double eta_jakobi = FastCaloSim::eta_jakobi(center_eta);

  const double delta_eta = delta_eta_mm / eta_jakobi / dist000;
  const double delta_phi = delta_phi_mm / center_r;

  hit.setEtaPhiZE(
      center_eta + delta_eta, center_phi + delta_phi, center_z, hit.E());

  FCS_MSG_DEBUG("HIT: E=" << hit.E() << " cs=" << cs << " eta=" << hit.eta()
                          << " phi=" << hit.phi() << " z=" << hit.z()
                          << " r=" << r << " alpha=" << alpha);

  return FCSSuccess;
}

auto TFCSHistoLateralShapeParametrization::Initialize(TH2* hist) -> bool
{
  if (hist == nullptr) {
    return false;
  }
  m_hist.Initialize(hist);
  if (m_hist.get_HistoContents().empty()) {
    return false;
  }

  set_number_of_hits(static_cast<float>(hist->Integral()));

  return true;
}

auto TFCSHistoLateralShapeParametrization::Initialize(const char* filepath,
                                                      const char* histname)
    -> bool
{
  // input file with histogram to fit
  std::unique_ptr<TFile> inputfile(TFile::Open(filepath, "READ"));
  if (inputfile == nullptr) {
    return false;
  }

  // histogram with hit pattern
  TH2* inputShape = static_cast<TH2*>(inputfile->Get(histname));
  if (inputShape == nullptr) {
    return false;
  }

  bool const OK = Initialize(inputShape);

  inputfile->Close();

  return OK;
}

void TFCSHistoLateralShapeParametrization::Print(Option_t* option) const
{
  const auto [shortprint, longprint, optprint] = print_options(option);
  TFCSLateralShapeParametrizationHitBase::Print(option);

  if (longprint) {
    if (is_phi_symmetric()) {
      FCS_MSG_INFO(optprint
                   << "  Histo: " << m_hist.get_HistoBordersx().size() - 1
                   << "*" << m_hist.get_HistoBordersy().size() - 1
                   << " bins, #hits=" << m_nhits << ", r scale=" << m_r_scale
                   << ", r offset=" << m_r_offset << "mm (phi symmetric)");
    } else {
      FCS_MSG_INFO(optprint
                   << "  Histo: " << m_hist.get_HistoBordersx().size() - 1
                   << "*" << m_hist.get_HistoBordersy().size() - 1
                   << " bins, #hits=" << m_nhits << ", r scale=" << m_r_scale
                   << ", r offset=" << m_r_offset << "mm (not phi symmetric)");
    }
  }
}
