// Copyright (c) 2024 CERN for the benefit of the FastCaloSim project

#include <iostream>
#include <limits>

#include "FastCaloSim/Core/TFCSEnergyAndHitGANV2.h"

#include <CLHEP/Random/RanluxEngine.h>

#include "CLHEP/Random/RandFlat.h"
#include "FastCaloSim/Core/TFCSExtrapolationState.h"
#include "FastCaloSim/Core/TFCSLateralShapeParametrizationHitBase.h"
#include "FastCaloSim/Core/TFCSSimulationState.h"
#include "FastCaloSim/Core/TFCSTruthState.h"
#include "FastCaloSim/Definitions/ParticleData.h"
#include "FastCaloSim/Geometry/CaloGeo.h"
#include "TF1.h"

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
  if (m_slice != nullptr) {
    delete m_slice;
  }
}

bool TFCSEnergyAndHitGANV2::is_match_calosample(int calosample) const
{
  if (get_Binning().find(calosample) == get_Binning().cend())
    return false;
  if (get_Binning().at(calosample).GetNbinsX() == 1)
    return false;
  return true;
}

unsigned int TFCSEnergyAndHitGANV2::get_nr_of_init(unsigned int bin) const
{
  if (bin >= m_bin_ninit.size())
    return 0;
  return m_bin_ninit[bin];
}

void TFCSEnergyAndHitGANV2::set_nr_of_init(unsigned int bin, unsigned int ninit)
{
  if (bin >= m_bin_ninit.size()) {
    m_bin_ninit.resize(bin + 1, 0);
    m_bin_ninit.shrink_to_fit();
  }
  m_bin_ninit[bin] = ninit;
}

// initialize lwtnn network
bool TFCSEnergyAndHitGANV2::initializeNetwork(
    int const& pid,
    int const& etaMin,
    const std::string& FastCaloGANInputFolderName)
{
  // initialize all necessary constants
  // FIXME eventually all these could be stored in the .json file

  FCS_MSG_INFO(
      "Using FastCaloGANInputFolderName: " << FastCaloGANInputFolderName);
  // get neural net JSON file as an std::istream object
  const int etaMax = etaMin + 5;

  reset_match_all_pdgid();
  set_pdgid(pid);
  if (pid == 11)
    add_pdgid(-pid);
  if (pid == 211)
    add_pdgid(-pid);
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

const std::string TFCSEnergyAndHitGANV2::get_variable_text(
    TFCSSimulationState& simulstate,
    const TFCSTruthState*,
    const TFCSExtrapolationState*) const
{
  return std::string(
      Form("layer=%d", simulstate.getAuxInfo<int>("GANlayer"_FCShash)));
}

bool TFCSEnergyAndHitGANV2::fillEnergy(
    TFCSSimulationState& simulstate,
    const TFCSTruthState* truth,
    const TFCSExtrapolationState* extrapol) const
{
  if (!truth) {
    FCS_MSG_ERROR("Invalid truth pointer");
    return false;
  }

  if (!extrapol) {
    FCS_MSG_ERROR("Invalid extrapolation pointer");
    return false;
  }

  const int pdgId = truth->pdgid();
  const float charge = ParticleData::charge(pdgId);

  float Einit;
  const float Ekin = truth->Ekin();
  if (OnlyScaleEnergy())
    Einit = simulstate.E();
  else
    Einit = Ekin;

  FCS_MSG_VERBOSE("Momentum " << truth->P() << " pdgId " << truth->pdgid());
  // check that the network exists
  if (!m_slice->IsGanCorrectlyLoaded()) {
    FCS_MSG_WARNING("GAN not loaded correctly.");
    return false;
  }

  const TFCSGANEtaSlice::NetworkOutputs& outputs =
      m_slice->GetNetworkOutputs(truth, extrapol, simulstate);
  FCS_MSG_VERBOSE("network outputs size: " << outputs.size());

  const TFCSGANXMLParameters::Binning& binsInLayers = m_param.GetBinning();
  const auto ganVersion = m_param.GetGANVersion();
  const TFCSGANEtaSlice::FitResultsPerLayer& fitResults =
      m_slice->GetFitResults();  // used only if GAN version > 1

  FCS_MSG_DEBUG("energy voxels size = " << outputs.size());

  double totalEnergy = 0;
  for (auto output : outputs) {
    totalEnergy += output.second;
  }
  if (totalEnergy < 0) {
    FCS_MSG_WARNING("Energy from GAN is negative, skipping particle");
    return false;
  }

  FCS_MSG_VERBOSE("Get binning");

  simulstate.set_E(0);

  int vox = 0;
  for (const auto& element : binsInLayers) {
    const int layer = element.first;
    const TH2D* h = &element.second;

    const int xBinNum = h->GetNbinsX();
    const int yBinNum = h->GetNbinsY();
    const TAxis* x = h->GetXaxis();

    // If only one bin in r means layer is empty, no value should be added
    if (xBinNum == 1) {
      FCS_MSG_VERBOSE(" Layer "
                      << layer
                      << " has only one bin in r, this means is it not used, "
                         "skipping (this is needed to keep correct "
                         "synchronization of voxel and layers)");
      // delete h;
      continue;
    }

    FCS_MSG_VERBOSE(" Getting energy for Layer " << layer);

    // First fill energies
    for (int ix = 1; ix <= xBinNum; ++ix) {
      double binsInAlphaInRBin = GetAlphaBinsForRBin(x, ix, yBinNum);
      for (int iy = 1; iy <= binsInAlphaInRBin; ++iy) {
        const double energyInVoxel = outputs.at(std::to_string(vox));
        FCS_MSG_VERBOSE(" Vox " << vox << " energy " << energyInVoxel
                                << " binx " << ix << " biny " << iy);

        if (energyInVoxel <= 0) {
          vox++;
          continue;
        }

        simulstate.add_E(layer, Einit * energyInVoxel);
        vox++;
      }
    }
  }

  for (unsigned int ichain = m_bin_start.back(); ichain < size(); ++ichain) {
    FCS_MSG_DEBUG("now run for all bins: " << chain()[ichain]->GetName());
    if (simulate_and_retry(chain()[ichain], simulstate, truth, extrapol)
        != FCSSuccess)
    {
      return FCSFatal;
    }
  }

  vox = 0;
  for (const auto& element : binsInLayers) {
    const int layer = element.first;
    const TH2D* h = &element.second;
    const int xBinNum = h->GetNbinsX();
    const int yBinNum = h->GetNbinsY();
    const TAxis* x = h->GetXaxis();
    const TAxis* y = h->GetYaxis();

    simulstate.setAuxInfo<int>("GANlayer"_FCShash, layer);
    TFCSLateralShapeParametrizationHitBase::Hit hit;

    // If only one bin in r means layer is empty, no value should be added
    if (xBinNum == 1) {
      FCS_MSG_VERBOSE(" Layer "
                      << layer
                      << " has only one bin in r, this means is it not used, "
                         "skipping (this is needed to keep correct "
                         "synchronization of voxel and layers)");
      // delete h;
      continue;
    }

    if (get_number_of_bins() > 0) {
      const int bin = get_bin(simulstate, truth, extrapol);
      if (bin >= 0 && bin < (int)get_number_of_bins()) {
        for (unsigned int ichain = m_bin_start[bin];
             ichain < TMath::Min(m_bin_start[bin] + get_nr_of_init(bin),
                                 m_bin_start[bin + 1]);
             ++ichain)
        {
          FCS_MSG_DEBUG("for " << get_variable_text(simulstate, truth, extrapol)
                               << " run init " << get_bin_text(bin) << ": "
                               << chain()[ichain]->GetName());
          if (chain()[ichain]->InheritsFrom(
                  TFCSLateralShapeParametrizationHitBase::Class()))
          {
            TFCSLateralShapeParametrizationHitBase* sim =
                (TFCSLateralShapeParametrizationHitBase*)(chain()[ichain]);
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
                        << ": " << get_bin_text(bin));
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

    const float dist000 =
        TMath::Sqrt(center_r * center_r + center_z * center_z);
    const float eta_jakobi = TMath::Abs(2.0 * TMath::Exp(-center_eta)
                                        / (1.0 + TMath::Exp(-2 * center_eta)));

    int nHitsAlpha;
    int nHitsR;

    // Now create hits
    for (int ix = 1; ix <= xBinNum; ++ix) {
      const int binsInAlphaInRBin = GetAlphaBinsForRBin(x, ix, yBinNum);

      // Horrible work around for variable # of bins along alpha direction
      const int binsToMerge = yBinNum == 32 ? 32 / binsInAlphaInRBin : 1;
      for (int iy = 1; iy <= binsInAlphaInRBin; ++iy) {
        const double energyInVoxel = outputs.at(std::to_string(vox));
        const int lowEdgeIndex = (iy - 1) * binsToMerge + 1;

        FCS_MSG_VERBOSE(" Vox " << vox << " energy " << energyInVoxel
                                << " binx " << ix << " biny " << iy);

        if (energyInVoxel <= 0) {
          vox++;
          continue;
        }

        if (fabs(pdgId) == 22 || fabs(pdgId) == 11) {
          // maximum 10 MeV per hit, equally distributed in alpha and r
          int maxHitsInVoxel = energyInVoxel * truth->Ekin() / 10;
          if (maxHitsInVoxel < 1)
            maxHitsInVoxel = 1;
          nHitsAlpha = sqrt(maxHitsInVoxel);
          nHitsR = sqrt(maxHitsInVoxel);
        } else {
          // One hit per mm along r
          nHitsR = x->GetBinUpEdge(ix) - x->GetBinLowEdge(ix);
          if (yBinNum == 1) {
            // nbins in alpha depend on circumference length
            const double r = x->GetBinUpEdge(ix);
            nHitsAlpha = ceil(2 * TMath::Pi() * r / binResolution);
          } else {
            // d = 2*r*sin (a/2r) this distance at the upper r must be 1mm for
            // layer 1 or 5, 5mm otherwise.
            const TAxis* y = h->GetYaxis();
            const double angle = y->GetBinUpEdge(iy) - y->GetBinLowEdge(iy);
            const double r = x->GetBinUpEdge(ix);
            const double d = 2 * r * sin(angle / 2 * r);
            nHitsAlpha = ceil(d / binResolution);
          }

          if (layer != 1 && layer != 5) {
            // For layers that are not EMB1 or EMEC1 use a maximum of 10 hits
            // per direction, a higher granularity is needed for the other
            // layers
            const int maxNhits = 10;
            nHitsAlpha = std::min(maxNhits, std::max(1, nHitsAlpha));
            nHitsR = std::min(maxNhits, std::max(1, nHitsR));
          }
        }

        for (int ir = 0; ir < nHitsR; ++ir) {
          double r =
              x->GetBinLowEdge(ix) + x->GetBinWidth(ix) / (nHitsR + 1) * ir;

          for (int ialpha = 1; ialpha <= nHitsAlpha; ++ialpha) {
            if (ganVersion > 1) {
              if (fitResults.at(layer)[ix - 1] != 0) {
                int tries = 0;
                double a = CLHEP::RandFlat::shoot(simulstate.randomEngine(),
                                                  x->GetBinLowEdge(ix),
                                                  x->GetBinUpEdge(ix));
                double rand_r =
                    log((a - x->GetBinLowEdge(ix)) / (x->GetBinWidth(ix)))
                    / fitResults.at(layer)[ix - 1];
                while ((rand_r < x->GetBinLowEdge(ix)
                        || rand_r > x->GetBinUpEdge(ix))
                       && tries < 100)
                {
                  a = CLHEP::RandFlat::shoot(simulstate.randomEngine(),
                                             x->GetBinLowEdge(ix),
                                             x->GetBinUpEdge(ix));
                  rand_r =
                      log((a - x->GetBinLowEdge(ix)) / (x->GetBinWidth(ix)))
                      / fitResults.at(layer)[ix - 1];
                  tries++;
                }
                if (tries >= 100) {
                  FCS_MSG_VERBOSE(" Too many tries for bin ["
                                  << x->GetBinLowEdge(ix) << "-"
                                  << x->GetBinUpEdge(ix) << "] having slope "
                                  << fitResults.at(layer)[ix - 1]
                                  << " will use grid (old method)");
                } else {
                  r = rand_r;
                }
              }
            }

            double alpha;
            if (binsInAlphaInRBin == 1) {
              alpha = CLHEP::RandFlat::shoot(
                  simulstate.randomEngine(), -TMath::Pi(), TMath::Pi());
            } else {
              alpha = y->GetBinLowEdge(lowEdgeIndex)
                  + y->GetBinWidth(iy) * binsToMerge / (nHitsAlpha + 1)
                      * ialpha;

              if (m_param.IsSymmetrisedAlpha()) {
                if (CLHEP::RandFlat::shoot(
                        simulstate.randomEngine(), -TMath::Pi(), TMath::Pi())
                    < 0)
                {
                  alpha = -alpha;
                }
              }
            }

            hit.reset();
            hit.set_E(Einit * energyInVoxel / (nHitsAlpha * nHitsR));

            if (layer <= 20) {
              float delta_eta_mm = r * cos(alpha);
              float delta_phi_mm = r * sin(alpha);

              FCS_MSG_VERBOSE("delta_eta_mm " << delta_eta_mm
                                              << " delta_phi_mm "
                                              << delta_phi_mm);

              // Particles with negative eta are expected to have the same shape
              // as those with positive eta after transformation: delta_eta -->
              // -delta_eta
              if (center_eta < 0.)
                delta_eta_mm = -delta_eta_mm;
              // We derive the shower shapes for electrons and positively
              // charged hadrons. Particle with the opposite charge are expected
              // to have the same shower shape after the transformation:
              // delta_phi --> -delta_phi
              if ((charge < 0. && pdgId != 11) || pdgId == -11)
                delta_phi_mm = -delta_phi_mm;

              const float delta_eta = delta_eta_mm / eta_jakobi / dist000;
              const float delta_phi = delta_phi_mm / center_r;

              hit.set_eta_x(center_eta + delta_eta);
              hit.set_phi_y(TVector2::Phi_mpi_pi(center_phi + delta_phi));

              FCS_MSG_VERBOSE(" Hit eta " << hit.eta() << " phi " << hit.phi()
                                          << " layer " << layer);
            } else {  // FCAL is in (x,y,z)
              const float hit_r = r * cos(alpha) + center_r;
              float delta_phi = r * sin(alpha) / center_r;
              // We derive the shower shapes for electrons and positively
              // charged hadrons. Particle with the opposite charge are expected
              // to have the same shower shape after the transformation:
              // delta_phi --> -delta_phi
              if ((charge < 0. && pdgId != 11) || pdgId == -11)
                delta_phi = -delta_phi;
              const float hit_phi =
                  TVector2::Phi_mpi_pi(center_phi + delta_phi);
              hit.set_eta_x(hit_r * cos(hit_phi));
              hit.set_phi_y(hit_r * sin(hit_phi));
              hit.set_z(center_z);
              FCS_MSG_VERBOSE(" Hit x " << hit.x() << " y " << hit.y()
                                        << " layer " << layer);
            }

            if (get_number_of_bins() > 0) {
              const int bin = get_bin(simulstate, truth, extrapol);
              if (bin >= 0 && bin < (int)get_number_of_bins()) {
                for (unsigned int ichain =
                         m_bin_start[bin] + get_nr_of_init(bin);
                     ichain < m_bin_start[bin + 1];
                     ++ichain)
                {
                  FCS_MSG_DEBUG(
                      "for " << get_variable_text(simulstate, truth, extrapol)
                             << " run " << get_bin_text(bin) << ": "
                             << chain()[ichain]->GetName());
                  if (chain()[ichain]->InheritsFrom(
                          TFCSLateralShapeParametrizationHitBase::Class()))
                  {
                    TFCSLateralShapeParametrizationHitBase* sim =
                        (TFCSLateralShapeParametrizationHitBase*)(chain()
                                                                      [ichain]);
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
                    << get_bin_text(bin));
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
    for (int ilayer = 0; ilayer < m_geo->n_layers(); ++ilayer) {
      simulstate.set_Efrac(ilayer, simulstate.E(ilayer) / simulstate.E());
    }
  }

  FCS_MSG_VERBOSE("Done particle");
  return true;
}

FCSReturnCode TFCSEnergyAndHitGANV2::simulate(
    TFCSSimulationState& simulstate,
    const TFCSTruthState* truth,
    const TFCSExtrapolationState* extrapol) const
{
  for (unsigned int ichain = 0; ichain < m_bin_start[0]; ++ichain) {
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

void TFCSEnergyAndHitGANV2::Print(Option_t* option) const
{
  TFCSParametrization::Print(option);
  TString opt(option);
  const bool shortprint = opt.Index("short") >= 0;
  const bool longprint =
      msgLvl(FCS_MSG::DEBUG) || (msgLvl(FCS_MSG::INFO) && !shortprint);
  TString optprint = opt;
  optprint.ReplaceAll("short", "");

  TString prefix = "- ";
  for (unsigned int ichain = 0; ichain < size(); ++ichain) {
    if (ichain == 0 && ichain != m_bin_start.front()) {
      prefix = "> ";
      if (longprint)
        FCS_MSG_INFO(optprint << prefix << "Run for all bins");
    }
    for (unsigned int ibin = 0; ibin < get_number_of_bins(); ++ibin) {
      if (ichain == m_bin_start[ibin]) {
        if (ibin < get_number_of_bins() - 1)
          if (ichain == m_bin_start[ibin + 1])
            continue;
        prefix = Form("%-2d", ibin);
        if (longprint)
          FCS_MSG_INFO(optprint << prefix << "Run for " << get_bin_text(ibin));
      }
    }
    if (ichain == m_bin_start.back()) {
      prefix = "< ";
      if (longprint)
        FCS_MSG_INFO(optprint << prefix << "Run for all bins");
    }
    chain()[ichain]->Print(opt + prefix);
  }
}

int TFCSEnergyAndHitGANV2::GetBinsInFours(double const& bins)
{
  if (bins < 4)
    return 4;
  else if (bins < 8)
    return 8;
  else if (bins < 16)
    return 16;
  else
    return 32;
}

int TFCSEnergyAndHitGANV2::GetAlphaBinsForRBin(const TAxis* x,
                                               int ix,
                                               int yBinNum) const
{
  double binsInAlphaInRBin = yBinNum;
  if (yBinNum == 32) {
    FCS_MSG_DEBUG("yBinNum is special value 32");
    const double widthX = x->GetBinWidth(ix);
    const double radious = x->GetBinCenter(ix);
    double circumference = radious * 2 * TMath::Pi();
    if (m_param.IsSymmetrisedAlpha()) {
      circumference = radious * TMath::Pi();
    }

    const double bins = circumference / widthX;
    binsInAlphaInRBin = GetBinsInFours(bins);
    FCS_MSG_DEBUG("Bin in alpha: " << binsInAlphaInRBin << " for r bin: " << ix
                                   << " (" << x->GetBinLowEdge(ix) << "-"
                                   << x->GetBinUpEdge(ix) << ")");
  }
  return binsInAlphaInRBin;
}
