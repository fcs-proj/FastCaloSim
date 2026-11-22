// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <format>
#include <iostream>
#include <limits>
#include <numbers>
#include <optional>
#include <vector>

#include "FastCaloSim/Core/TFCSEnergyAndHitGANV2.h"

#include <CLHEP/Random/RanluxEngine.h>

#include "CLHEP/Random/RandFlat.h"
#include "FastCaloSim/Core/TFCSExtrapolationState.h"
#include "FastCaloSim/Core/TFCSLateralShapeParametrizationHitBase.h"
#include "FastCaloSim/Core/TFCSSimulationState.h"
#include "FastCaloSim/Core/TFCSTruthState.h"
#include "FastCaloSim/Definitions/FloatCompare.h"
#include "FastCaloSim/Definitions/ParticleData.h"
#include "FastCaloSim/Geometry/CaloGeo.h"
#include "TF1.h"

namespace
{
// Draws the radius of a hit in the r bin [low, low + width] from the
// exponential exp(slope * r) fitted to the radial profile of the bin,
// truncated to the bin, with one uniform random number by inverting its
// distribution function. With k = |slope| and x the distance from the edge
// towards which the exponential falls (the lower edge for a negative slope,
// the upper edge for a positive one), F(x) = (1 - exp(-k x)) /
// (1 - exp(-k width)), so x = -log1p(v * expm1(-k width)) / k. In this form
// expm1 and log1p stay accurate for small k * width (where the distribution
// becomes uniform) and nothing overflows for large ones. A slope of zero (no
// fit, or GAN version 1) or a non-finite one describes no shape: such bins
// keep the grid radius.
class RadialSampler
{
public:
  RadialSampler(double slope, double low, double width)
      : m_slope(slope)
      , m_edge(slope < 0 ? low : low + width)
      , m_expm1_width(std::expm1(-std::abs(slope) * width))
      , m_sampled(std::isfinite(slope)
                  && !FastCaloSim::exactly_equal(slope, 0.0))
  {
  }

  // The sampled radius, or nothing to keep the grid radius
  auto sample(CLHEP::HepRandomEngine* engine) const -> std::optional<double>
  {
    if (!m_sampled) {
      return std::nullopt;
    }
    const double v = CLHEP::RandFlat::shoot(engine);
    // log1p(...) / m_slope is the signed step from m_edge into the bin
    return m_edge + (std::log1p(v * m_expm1_width) / m_slope);
  }

private:
  double m_slope;
  double m_edge;
  double m_expm1_width;
  bool m_sampled;
};
}  // namespace

//=============================================
//======= TFCSEnergyAndHitGANV2 =========
//=============================================

TFCSEnergyAndHitGANV2::TFCSEnergyAndHitGANV2(const char* name,
                                             const char* title,
                                             CaloGeo* geo)
    : TFCSParametrizationBinnedChain(name, title)
    , m_geo(geo)
{
  set_GANfreemem();
}

TFCSEnergyAndHitGANV2::~TFCSEnergyAndHitGANV2()
{
  delete m_slice;
}

auto TFCSEnergyAndHitGANV2::is_match_calosample(std::size_t calosample) const
    -> bool
{
  // The binning is stored per layer with int keys
  const auto it = get_Binning().find(static_cast<int>(calosample));
  return it != get_Binning().cend() && it->second.GetNbinsX() != 1;
}

auto TFCSEnergyAndHitGANV2::get_nr_of_init(std::size_t bin) const
    -> unsigned int
{
  if (bin >= m_bin_ninit.size()) {
    return 0;
  }
  return static_cast<unsigned int>(m_bin_ninit[bin]);
}

void TFCSEnergyAndHitGANV2::set_nr_of_init(std::size_t bin, unsigned int ninit)
{
  if (bin >= m_bin_ninit.size()) {
    m_bin_ninit.resize(bin + 1, 0);
    m_bin_ninit.shrink_to_fit();
  }
  m_bin_ninit[bin] = static_cast<int>(ninit);
}

// initialize lwtnn network
auto TFCSEnergyAndHitGANV2::initializeNetwork(
    int const& pid,
    int const& etaMin,
    const std::string& FastCaloGANInputFolderName) -> bool
{
  // initialize all necessary constants
  // FIXME eventually all these could be stored in the .json file

  FCS_MSG_INFO(
      "Using FastCaloGANInputFolderName: " << FastCaloGANInputFolderName);
  // get neural net JSON file as an std::istream object
  const int etaMax = etaMin + 5;

  reset_match_all_pdgid();
  set_pdgid(pid);
  if (pid == 11) {
    add_pdgid(-pid);
  }
  if (pid == 211) {
    add_pdgid(-pid);
  }
  set_eta_min(etaMin / 100.0);
  set_eta_max(etaMax / 100.0);
  set_eta_nominal((etaMin + etaMax) / 200.0);

  int pidForXml = pid;
  if (pid != 22 && pid != 11) {
    pidForXml = 211;
  }

  const int etaMid = (etaMin + etaMax) / 2;
  m_param.InitialiseFromXML(pidForXml, etaMid, FastCaloGANInputFolderName);
  m_param.Print();
  m_slice = new TFCSGANEtaSlice(pid, etaMin, etaMax, m_param);
  m_slice->Print();
  return m_slice->LoadGAN();
}

auto TFCSEnergyAndHitGANV2::get_variable_text(
    TFCSSimulationState& simulstate,
    const TFCSTruthState* /*unused*/,
    const TFCSExtrapolationState* /*unused*/) const -> const std::string
{
  return std::format("layer={}",
                     simulstate.getAuxInfo<int>("GANlayer"_FCShash));
}

auto TFCSEnergyAndHitGANV2::fillEnergy(
    TFCSSimulationState& simulstate,
    const TFCSTruthState* truth,
    const TFCSExtrapolationState* extrapol) const -> bool
{
  if (truth == nullptr) {
    FCS_MSG_ERROR("Invalid truth pointer");
    return false;
  }

  if (extrapol == nullptr) {
    FCS_MSG_ERROR("Invalid extrapolation pointer");
    return false;
  }

  const int pdgId = truth->pdgid();
  const double charge = ParticleData::charge(pdgId);

  const double Ekin = truth->Ekin();
  const double Einit = OnlyScaleEnergy() ? simulstate.E() : Ekin;

  FCS_MSG_VERBOSE("Momentum " << truth->P() << " pdgId " << truth->pdgid());
  // check that the network exists
  if (!m_slice->IsGanCorrectlyLoaded()) {
    FCS_MSG_WARNING("GAN not loaded correctly.");
    return false;
  }

  const std::vector<double> voxelFractions =
      m_slice->GetVoxelEnergyFractions(truth, extrapol, simulstate);
  if (voxelFractions.empty()) {
    FCS_MSG_ERROR("No voxel energies from the GAN");
    return false;
  }

  const TFCSGANXMLParameters::Binning& binsInLayers = m_param.GetBinning();
  const auto ganVersion = m_param.GetGANVersion();
  const TFCSGANEtaSlice::FitResultsPerLayer& fitResults =
      m_slice->GetFitResults();  // used only if GAN version > 1

  FCS_MSG_DEBUG("energy voxels size = " << voxelFractions.size());

  double totalEnergy = 0;
  for (const double fraction : voxelFractions) {
    totalEnergy += fraction;
  }
  if (totalEnergy < 0) {
    FCS_MSG_WARNING("Energy from GAN is negative, skipping particle");
    return false;
  }

  FCS_MSG_VERBOSE("Get binning");

  simulstate.set_E(0);

  std::size_t vox = 0;
  for (const auto& [layer, h] : binsInLayers) {
    const int xBinNum = h.GetNbinsX();
    const int yBinNum = h.GetNbinsY();
    const TAxis* x = h.GetXaxis();

    // If only one bin in r means layer is empty, no value should be added
    if (xBinNum == 1) {
      FCS_MSG_VERBOSE(" Layer "
                      << layer
                      << " has only one bin in r, this means is it not used, "
                         "skipping (this is needed to keep correct "
                         "synchronization of voxel and layers)");
      continue;
    }

    FCS_MSG_VERBOSE(" Getting energy for Layer " << layer);

    // First fill energies
    for (int ix = 1; ix <= xBinNum; ++ix) {
      double const binsInAlphaInRBin = GetAlphaBinsForRBin(x, ix, yBinNum);
      for (int iy = 1; iy <= binsInAlphaInRBin; ++iy) {
        const double energyInVoxel = voxelFractions.at(vox);
        FCS_MSG_VERBOSE(" Vox " << vox << " energy " << energyInVoxel
                                << " binx " << ix << " biny " << iy);

        if (energyInVoxel <= 0) {
          vox++;
          continue;
        }

        simulstate.add_E(static_cast<std::size_t>(layer),
                         Einit * energyInVoxel);
        vox++;
      }
    }
  }

  for (std::size_t ichain = m_bin_start.back(); ichain < size(); ++ichain) {
    FCS_MSG_DEBUG("now run for all bins: " << chain()[ichain]->GetName());
    if (simulate_and_retry(chain()[ichain], simulstate, truth, extrapol)
        != FCSSuccess)
    {
      return false;
    }
  }

  vox = 0;
  for (const auto& [layer, h] : binsInLayers) {
    const int xBinNum = h.GetNbinsX();
    const int yBinNum = h.GetNbinsY();
    const TAxis* x = h.GetXaxis();
    const TAxis* y = h.GetYaxis();

    simulstate.setAuxInfo<int>("GANlayer"_FCShash, layer);
    TFCSLateralShapeParametrizationHitBase::Hit hit;

    // If only one bin in r means layer is empty, no value should be added
    if (xBinNum == 1) {
      FCS_MSG_VERBOSE(" Layer "
                      << layer
                      << " has only one bin in r, this means is it not used, "
                         "skipping (this is needed to keep correct "
                         "synchronization of voxel and layers)");
      continue;
    }

    if (get_number_of_bins() > 0) {
      // Select once: the selection may draw random numbers
      const auto selected = get_bin(simulstate, truth, extrapol);
      const auto bin = valid_bin(selected);
      if (bin) {
        const std::size_t ubin = *bin;
        for (unsigned int ichain = m_bin_start[ubin];
             ichain < std::min(m_bin_start[ubin] + get_nr_of_init(ubin),
                               m_bin_start[ubin + 1]);
             ++ichain)
        {
          FCS_MSG_DEBUG("for " << get_variable_text(simulstate, truth, extrapol)
                               << " run init " << get_bin_text(bin) << ": "
                               << chain()[ichain]->GetName());
          if (chain()[ichain]->InheritsFrom(
                  TFCSLateralShapeParametrizationHitBase::Class()))
          {
            auto* sim = static_cast<TFCSLateralShapeParametrizationHitBase*>(
                chain()[ichain]);
            if (sim->simulate_hit(hit, simulstate, truth, extrapol)
                != FCSSuccess)
            {
              FCS_MSG_ERROR("error for "
                            << get_variable_text(simulstate, truth, extrapol)
                            << " run init " << get_bin_text(bin) << ": "
                            << chain()[ichain]->GetName());
              return false;
            }
          } else {
            FCS_MSG_ERROR("for "
                          << get_variable_text(simulstate, truth, extrapol)
                          << " run init " << get_bin_text(bin) << ": "
                          << chain()[ichain]->GetName()
                          << " does not inherit from "
                             "TFCSLateralShapeParametrizationHitBase");
            return false;
          }
        }
      } else {
        FCS_MSG_WARNING("nothing to init for "
                        << get_variable_text(simulstate, truth, extrapol)
                        << ": " << get_bin_text(selected));
      }
    }

    int binResolution = 5;
    if (layer == 1 || layer == 5) {
      binResolution = 1;
    }

    const double center_eta = hit.center_eta();
    const double center_phi = hit.center_phi();
    const double center_r = hit.center_r();
    const double center_z = hit.center_z();

    FCS_MSG_VERBOSE(" Layer " << layer << " Extrap eta " << center_eta
                              << " phi " << center_phi << " R " << center_r);

    const double dist000 = FastCaloSim::dist000(center_r, center_z);
    const double eta_jakobi = FastCaloSim::eta_jakobi(center_eta);

    int nHitsAlpha {};
    int nHitsR {};

    // Now create hits
    for (int ix = 1; ix <= xBinNum; ++ix) {
      const int binsInAlphaInRBin = GetAlphaBinsForRBin(x, ix, yBinNum);
      const auto r_bin = static_cast<std::size_t>(ix - 1);
      // From GAN version 2 on, the radii are sampled within the r bin
      const RadialSampler radial(
          ganVersion > 1 ? fitResults.at(layer)[r_bin] : 0.0,
          x->GetBinLowEdge(ix),
          x->GetBinWidth(ix));

      // Horrible work around for variable # of bins along alpha direction
      const int binsToMerge = yBinNum == 32 ? 32 / binsInAlphaInRBin : 1;
      for (int iy = 1; iy <= binsInAlphaInRBin; ++iy) {
        const double energyInVoxel = voxelFractions.at(vox);
        const int lowEdgeIndex = (iy - 1) * binsToMerge + 1;
        // The voxel spans binsToMerge equal alpha bins (radians)
        const double alphaWidth = y->GetBinWidth(iy) * binsToMerge;

        FCS_MSG_VERBOSE(" Vox " << vox << " energy " << energyInVoxel
                                << " binx " << ix << " biny " << iy);

        if (energyInVoxel <= 0) {
          vox++;
          continue;
        }

        if (std::abs(pdgId) == 22 || std::abs(pdgId) == 11) {
          // maximum 10 MeV per hit, equally distributed in alpha and r
          int maxHitsInVoxel =
              static_cast<int>(energyInVoxel * truth->Ekin() / 10);
          maxHitsInVoxel = std::max(maxHitsInVoxel, 1);
          nHitsAlpha = static_cast<int>(std::sqrt(maxHitsInVoxel));
          nHitsR = static_cast<int>(std::sqrt(maxHitsInVoxel));
        } else {
          // One hit per mm along r
          nHitsR = static_cast<int>(x->GetBinUpEdge(ix) - x->GetBinLowEdge(ix));
          if (yBinNum == 1) {
            // nbins in alpha depend on circumference length
            const double r = x->GetBinUpEdge(ix);
            nHitsAlpha = static_cast<int>(
                std::ceil(2 * std::numbers::pi * r / binResolution));
          } else {
            // One hit per 1 mm (layers 1 and 5) or 5 mm of the chord
            // d = 2 r sin(alphaWidth / 2) that the voxel spans at its upper r
            const double r = x->GetBinUpEdge(ix);
            const double d = 2 * r * std::sin(alphaWidth / 2);
            nHitsAlpha = static_cast<int>(std::ceil(d / binResolution));
          }

          // At least one hit in each direction, so that the energy of the
          // voxel, which is already part of the layer energy, is deposited
          nHitsAlpha = std::max(1, nHitsAlpha);
          nHitsR = std::max(1, nHitsR);

          if (layer != 1 && layer != 5) {
            // For layers that are not EMB1 or EMEC1 use a maximum of 10 hits
            // per direction, a higher granularity is needed for the other
            // layers
            const int maxNhits = 10;
            nHitsAlpha = std::min(maxNhits, nHitsAlpha);
            nHitsR = std::min(maxNhits, nHitsR);
          }
        }

        for (int ir = 0; ir < nHitsR; ++ir) {
          // Interior points k / (nHitsR + 1), k = 1..nHitsR, of the r bin, like
          // the alpha grid
          const double r_grid = x->GetBinLowEdge(ix)
              + x->GetBinWidth(ix) / (nHitsR + 1) * (ir + 1);

          for (int ialpha = 1; ialpha <= nHitsAlpha; ++ialpha) {
            const double r =
                radial.sample(simulstate.randomEngine()).value_or(r_grid);

            const bool random_alpha = binsInAlphaInRBin == 1;
            double alpha = random_alpha
                ? CLHEP::RandFlat::shoot(simulstate.randomEngine(),
                                         -std::numbers::pi,
                                         std::numbers::pi)
                : y->GetBinLowEdge(lowEdgeIndex)
                    + alphaWidth / (nHitsAlpha + 1) * ialpha;
            if (!random_alpha && m_param.IsSymmetrisedAlpha()) {
              if (CLHEP::RandFlat::shoot(simulstate.randomEngine(),
                                         -std::numbers::pi,
                                         std::numbers::pi)
                  < 0)
              {
                alpha = -alpha;
              }
            }

            hit.reset();
            hit.set_E(Einit * energyInVoxel / (nHitsAlpha * nHitsR));

            if (layer <= 20) {
              double delta_eta_mm = r * std::cos(alpha);
              double delta_phi_mm = r * std::sin(alpha);

              FCS_MSG_VERBOSE("delta_eta_mm " << delta_eta_mm
                                              << " delta_phi_mm "
                                              << delta_phi_mm);

              // Particles with negative eta are expected to have the same shape
              // as those with positive eta after transformation: delta_eta -->
              // -delta_eta
              if (center_eta < 0.) {
                delta_eta_mm = -delta_eta_mm;
              }
              // We derive the shower shapes for electrons and positively
              // charged hadrons. Particle with the opposite charge are expected
              // to have the same shower shape after the transformation:
              // delta_phi --> -delta_phi
              if ((charge < 0. && pdgId != 11) || pdgId == -11) {
                delta_phi_mm = -delta_phi_mm;
              }

              const double delta_eta = delta_eta_mm / eta_jakobi / dist000;
              const double delta_phi = delta_phi_mm / center_r;

              hit.set_eta_x(center_eta + delta_eta);
              hit.set_phi_y(TVector2::Phi_mpi_pi(center_phi + delta_phi));

              FCS_MSG_VERBOSE(" Hit eta " << hit.eta() << " phi " << hit.phi()
                                          << " layer " << layer);
            } else {  // FCAL is in (x,y,z)
              const double hit_r = r * std::cos(alpha) + center_r;
              double delta_phi = r * std::sin(alpha) / center_r;
              // We derive the shower shapes for electrons and positively
              // charged hadrons. Particle with the opposite charge are expected
              // to have the same shower shape after the transformation:
              // delta_phi --> -delta_phi
              if ((charge < 0. && pdgId != 11) || pdgId == -11) {
                delta_phi = -delta_phi;
              }
              const double hit_phi =
                  TVector2::Phi_mpi_pi(center_phi + delta_phi);
              hit.set_eta_x(hit_r * std::cos(hit_phi));
              hit.set_phi_y(hit_r * std::sin(hit_phi));
              hit.set_z(center_z);
              FCS_MSG_VERBOSE(" Hit x " << hit.x() << " y " << hit.y()
                                        << " layer " << layer);
            }

            if (get_number_of_bins() > 0) {
              // Select once: the selection may draw random numbers
              const auto selected = get_bin(simulstate, truth, extrapol);
              const auto bin = valid_bin(selected);
              if (bin) {
                const std::size_t ubin = *bin;
                for (unsigned int ichain =
                         m_bin_start[ubin] + get_nr_of_init(ubin);
                     ichain < m_bin_start[ubin + 1];
                     ++ichain)
                {
                  FCS_MSG_DEBUG(
                      "for " << get_variable_text(simulstate, truth, extrapol)
                             << " run " << get_bin_text(bin) << ": "
                             << chain()[ichain]->GetName());
                  if (chain()[ichain]->InheritsFrom(
                          TFCSLateralShapeParametrizationHitBase::Class()))
                  {
                    auto* sim =
                        static_cast<TFCSLateralShapeParametrizationHitBase*>(
                            chain()[ichain]);
                    if (sim->simulate_hit(hit, simulstate, truth, extrapol)
                        != FCSSuccess)
                    {
                      FCS_MSG_ERROR(
                          "error for "
                          << get_variable_text(simulstate, truth, extrapol)
                          << " run init " << get_bin_text(bin) << ": "
                          << chain()[ichain]->GetName());
                      return false;
                    }
                  } else {
                    FCS_MSG_ERROR(
                        "for " << get_variable_text(simulstate, truth, extrapol)
                               << " run init " << get_bin_text(bin) << ": "
                               << chain()[ichain]->GetName()
                               << " does not inherit from "
                                  "TFCSLateralShapeParametrizationHitBase");
                    return false;
                  }
                }
              } else {
                FCS_MSG_WARNING(
                    "nothing to do for "
                    << get_variable_text(simulstate, truth, extrapol) << ": "
                    << get_bin_text(selected));
              }
            } else {
              FCS_MSG_WARNING("no bins defined, is this intended?");
            }
          }
        }
        vox++;
      }
    }

    FCS_MSG_VERBOSE("Number of voxels " << vox);
    FCS_MSG_VERBOSE("Done layer " << layer);
  }

  if (simulstate.E() > std::numeric_limits<double>::epsilon()) {
    for (std::size_t ilayer = 0; ilayer < m_geo->n_layers(); ++ilayer) {
      simulstate.set_Efrac(ilayer, simulstate.E(ilayer) / simulstate.E());
    }
  }

  FCS_MSG_VERBOSE("Done particle");
  return true;
}

auto TFCSEnergyAndHitGANV2::simulate(
    TFCSSimulationState& simulstate,
    const TFCSTruthState* truth,
    const TFCSExtrapolationState* extrapol) const -> FCSReturnCode
{
  for (std::size_t ichain = 0; ichain < m_bin_start[0]; ++ichain) {
    FCS_MSG_DEBUG("now run for all bins: " << chain()[ichain]->GetName());
    if (simulate_and_retry(chain()[ichain], simulstate, truth, extrapol)
        != FCSSuccess)
    {
      return FCSFatal;
    }
  }

  FCS_MSG_VERBOSE("Fill Energies");
  if (!fillEnergy(simulstate, truth, extrapol)) {
    FCS_MSG_WARNING("Could not fill energies ");
    // bail out but do not stop the job
    return FCSSuccess;
  }

  return FCSSuccess;
}

auto TFCSEnergyAndHitGANV2::GetBinsInFours(double const bins) -> int
{
  if (bins < 4) {
    return 4;
  }
  if (bins < 8) {
    return 8;
  }
  if (bins < 16) {
    return 16;
  }
  return 32;
}

auto TFCSEnergyAndHitGANV2::GetAlphaBinsForRBin(const TAxis* x,
                                                int ix,
                                                int yBinNum) const -> int
{
  double binsInAlphaInRBin = yBinNum;
  if (yBinNum == 32) {
    FCS_MSG_DEBUG("yBinNum is special value 32");
    const double widthX = x->GetBinWidth(ix);
    const double radious = x->GetBinCenter(ix);
    double circumference = radious * 2. * std::numbers::pi;
    if (m_param.IsSymmetrisedAlpha()) {
      circumference = radious * std::numbers::pi;
    }

    const double bins = circumference / widthX;
    binsInAlphaInRBin = GetBinsInFours(bins);
    FCS_MSG_DEBUG("Bin in alpha: " << binsInAlphaInRBin << " for r bin: " << ix
                                   << " (" << x->GetBinLowEdge(ix) << "-"
                                   << x->GetBinUpEdge(ix) << ")");
  }
  return static_cast<int>(binsInAlphaInRBin);
}

void TFCSEnergyAndHitGANV2::fixHists()
{
  m_param.fixHists();
  // m_slice is persistified as well and holds its own copy of the binning
  // histograms, so it needs the same treatment. See ATLASSIM-7031.
  if (m_slice != nullptr) {
    m_slice->fixHists();
  }
}
