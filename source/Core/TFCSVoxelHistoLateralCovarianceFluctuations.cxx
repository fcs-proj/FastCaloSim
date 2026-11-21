// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#include <cmath>
#include <cstddef>
#include <format>
#include <memory>
#include <numbers>

#include "FastCaloSim/Core/TFCSVoxelHistoLateralCovarianceFluctuations.h"

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGauss.h"
#include "FastCaloSim/Core/TFCS1DFunction.h"
#include "FastCaloSim/Core/TFCS1DFunctionTemplateHistogram.h"
#include "FastCaloSim/Core/TFCSSimulationState.h"
#include "FastCaloSim/Geometry/CaloGeo.h"
#include "TFile.h"
#include "TH2.h"
#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"
#include "TVector2.h"

//=============================================
//======= TFCSVoxelHistoLateralCovarianceFluctuations =========
//=============================================

namespace
{
// The layer numbers are small, so they are valid AuxInfo keys
auto layer_hash(std::size_t layer) -> std::uint32_t
{
  return static_cast<std::uint32_t>(layer);
}
auto layer_hash_geo(std::size_t layer) -> std::uint32_t
{
  return static_cast<std::uint32_t>(layer) | 0x80000000u;
}
}  // namespace

TFCSVoxelHistoLateralCovarianceFluctuations::
    TFCSVoxelHistoLateralCovarianceFluctuations(const char* name,
                                                const char* title)
    : TFCSLateralShapeParametrizationHitBase(name, title)
{
}

TFCSVoxelHistoLateralCovarianceFluctuations::
    ~TFCSVoxelHistoLateralCovarianceFluctuations() = default;

auto TFCSVoxelHistoLateralCovarianceFluctuations::initialize(
    TFile* inputfile, const std::string& folder) -> bool
{
  // load m_eigenvariances and m_parMeans from input file
  // load histograms for each cell from input file
  // m_transform[i][j]=new TFCS1DFunctionInt16Int16Histogram(hist)
  // where hist is the 1D histogram from which a CDF was calculated

  inputfile->cd(folder.c_str());
  std::size_t cs = calosample();
  int bin = 1;

  TH2* temp = dynamic_cast<TH2*>(inputfile->Get(
      std::format("voxel_template_cs{}_pca{}", cs, bin).c_str()));
  if (!temp) {
    FCS_MSG_ERROR("Template hist not found for cs " + std::to_string(cs));
    return false;
  }
  m_voxel_template.push_back(temp);
  m_nDim_x = m_voxel_template[0]->GetNbinsX();
  m_nDim_y = m_voxel_template[0]->GetNbinsY();

  // TFile::Get returns a new object, owned by the caller, for objects that are
  // not histograms or trees
  while (true) {
    TVectorD parMeans;
    TMatrixD EigenVectors;
    TVectorD EigenValues;
    std::vector<std::vector<TFCS1DFunction*>> transform;

    std::string label = std::format("_cs{}_pca{}", cs, bin);

    std::unique_ptr<TObject> obj(inputfile->Get(("parMeans" + label).c_str()));
    if (!obj)
      break;
    if (msgLvl(FastCaloSim::MSG::DEBUG)) {
      FCS_MSG_DEBUG("parMeans");
      obj->Print();
    }
    parMeans.ResizeTo(m_nDim_x * m_nDim_y);
    parMeans = *dynamic_cast<TVectorT<double>*>(obj.get());
    m_parMeans.push_back(parMeans);

    obj.reset(inputfile->Get(("covMatrix" + label).c_str()));
    if (!obj) {
      FCS_MSG_ERROR("covMatrix" + label + " not found");
      return false;
    }
    TMatrixTSym<double> covMatrix =
        *dynamic_cast<TMatrixTSym<double>*>(obj.get());
    TMatrixDSymEigen eigenvariances = TMatrixDSymEigen(covMatrix);
    EigenVectors.ResizeTo(eigenvariances.GetEigenVectors());
    EigenVectors = eigenvariances.GetEigenVectors();
    m_EigenVectors.push_back(EigenVectors);

    EigenValues.ResizeTo(eigenvariances.GetEigenValues());
    EigenValues = eigenvariances.GetEigenValues();
    m_EigenValues.push_back(EigenValues);

    if (msgLvl(FastCaloSim::MSG::DEBUG)) {
      FCS_MSG_DEBUG("eigenvariances");
      eigenvariances.GetEigenValues().Print();
      eigenvariances.GetEigenVectors().Print();
    }
    if (msgLvl(FastCaloSim::MSG::DEBUG)) {
      FCS_MSG_DEBUG("m_EigenValues");
      EigenValues.Print();
      FCS_MSG_DEBUG("m_EigenVectors");
      EigenVectors.Print();
    }

    const auto n_x = static_cast<std::size_t>(m_nDim_x);
    const auto n_y = static_cast<std::size_t>(m_nDim_y);
    transform.resize(n_x);
    for (std::size_t x = 0; x < n_x; ++x) {
      transform[x].resize(n_y);
      for (std::size_t y = 0; y < n_y; ++y) {
        // Naming of hists is hist_xy where x, y are indices of cells
        std::string histname = std::format("hist_{}_{}{}", x, y, label);
        TH1* hist = static_cast<TH1*>(inputfile->Get(histname.c_str()));
        if (!hist) {
          FCS_MSG_ERROR("Histogram " << histname << " not found");
          return false;
        }
        TFCS1DFunctionInt32Int32Histogram* func =
            new TFCS1DFunctionInt32Int32Histogram();
        func->Initialize(hist, false);
        transform[x][y] = func;
      }
    }
    m_transform.push_back(std::move(transform));
    ++bin;
  }

  return true;
}

// Implementation of multivariate Gaussian with ROOT classes (from ROOT forum,
// but validated)
void TFCSVoxelHistoLateralCovarianceFluctuations::MultiGaus(
    TFCSSimulationState& simulstate, TVectorD& genPars) const
{
  const auto ibin = static_cast<std::size_t>(simulstate.Ebin() - 1);

  int nPars = m_parMeans[ibin].GetNrows();
  genPars.ResizeTo(nPars);

  TVectorD rotParMeans = m_EigenVectors[ibin] * m_parMeans[ibin];
  for (int iPar = 0; iPar < nPars; iPar++) {
    double variance = m_EigenValues[ibin][iPar];
    // check for positive-definiteness of covMatrix
    if (variance < 0) {
      FCS_MSG_ERROR("Got a negative eigenvariance (" << variance
                                                     << ") on iPar = " << iPar);
      variance = 0;
    }
    genPars[iPar] = CLHEP::RandGauss::shoot(
        simulstate.randomEngine(), rotParMeans[iPar], std::sqrt(variance));
    FCS_MSG_DEBUG("genPars[" << iPar << "]=" << genPars[iPar]
                             << " rotParMeans[iPar]=" << rotParMeans[iPar]
                             << " sqrt(variance)=" << std::sqrt(variance));
  }
  genPars = m_EigenVectors[ibin] * genPars;
}

auto TFCSVoxelHistoLateralCovarianceFluctuations::simulate(
    TFCSSimulationState& simulstate,
    const TFCSTruthState* /*truth*/,
    const TFCSExtrapolationState* /*extrapol*/) const -> FCSReturnCode
{
  if (!simulstate.randomEngine()) {
    return FCSFatal;
  }

  if (msgLvl(FastCaloSim::MSG::DEBUG)) {
    FCS_MSG_DEBUG("simulstate before clearing AuxInfo");
    simulstate.Print();
  }

  weight_t* weightvec;
  TH2* voxel_temp;

  int Ebin = simulstate.Ebin();

  for (std::size_t ilayer = 0; ilayer < m_geo->n_layers(); ++ilayer) {
    if (simulstate.hasAuxInfo(layer_hash(ilayer))) {
      weightvec = static_cast<weight_t*>(
          simulstate.getAuxInfo<void*>(layer_hash(ilayer)));
      if (weightvec) {
        delete weightvec;
        simulstate.setAuxInfo<void*>(layer_hash(ilayer), nullptr);
      }
    }
    if (simulstate.hasAuxInfo(layer_hash_geo(ilayer))) {
      voxel_temp = static_cast<TH2*>(
          simulstate.getAuxInfo<void*>(layer_hash_geo(ilayer)));
      if (voxel_temp) {
        delete voxel_temp;
        simulstate.setAuxInfo<void*>(layer_hash_geo(ilayer), nullptr);
      }
    }
  }

  if (msgLvl(FastCaloSim::MSG::DEBUG)) {
    FCS_MSG_DEBUG("simulstate after clearing AuxInfo");
    simulstate.Print();
  }
  if (m_parMeans.empty() || Ebin <= 0) {
    // not initialized, do nothing in case this instance is just used to clean
    // up memory
    return FCSSuccess;
  }

  if (static_cast<size_t>(Ebin) > m_parMeans.size()) {
    FCS_MSG_ERROR("Ebin " << Ebin << " exceeds loaded PCA tables (size "
                          << m_parMeans.size() << ")");
    return FCSFatal;
  }

  // TODO: the following code should be executed for all relevant calo layers,
  // possibly simulating correlated fluctuations between layers depending on the
  // PCA bin
  const auto ibin = static_cast<std::size_t>(Ebin - 1);
  const auto n_x = static_cast<std::size_t>(m_nDim_x);
  const auto n_y = static_cast<std::size_t>(m_nDim_y);
  weightvec = new weight_t(n_x);
  for (std::size_t x = 0; x < n_x; ++x)
    (*weightvec)[x].resize(n_y);

  TVectorD genPars(m_parMeans[ibin].GetNrows());

  // Fill genPars with multivarate random Gauss
  MultiGaus(simulstate, genPars);

  // Loop through voxels
  int count = 0;
  for (std::size_t x = 0; x < n_x; ++x) {
    for (std::size_t y = 0; y < n_y; ++y) {
      // Get value of CDF corresponding to generated Gauss random
      double cdf_val =
          (std::erf(genPars[count] * 2.0 / std::numbers::pi) + 1) / 2.0;

      // Map that cdf value to input bin and value
      double orig_val = m_transform[ibin][x][y]->rnd_to_fct(cdf_val);

      // Fill hist, keep linear count
      (*weightvec)[x][y] = orig_val;
      ++count;

      FCS_MSG_DEBUG("CELL x=" << x << " y=" << y << " : cdf_val=" << cdf_val
                              << " orig_val=" << orig_val);
    }
  }

  voxel_temp = static_cast<TH2*>(m_voxel_template[0]->Clone());

  // For now simulating only the layer calosample()
  const std::size_t cs = calosample();
  simulstate.setAuxInfo<void*>(layer_hash(cs), weightvec);
  simulstate.setAuxInfo<void*>(layer_hash_geo(cs), voxel_temp);

  if (msgLvl(FastCaloSim::MSG::DEBUG)) {
    FCS_MSG_DEBUG("simulstate after storing weight " << weightvec
                                                     << " in AuxInfo");
    simulstate.Print();
  }

  return FCSSuccess;
}

auto TFCSVoxelHistoLateralCovarianceFluctuations::simulate_hit(
    Hit& hit,
    TFCSSimulationState& simulstate,
    const TFCSTruthState* /*truth*/,
    const TFCSExtrapolationState* /*extrapol*/) -> FCSReturnCode
{
  const double center_eta = hit.center_eta();
  const double center_phi = hit.center_phi();
  const double center_r = hit.center_r();
  const double center_z = hit.center_z();

  const std::size_t cs = calosample();
  weight_t* weightvec = nullptr;
  if (simulstate.hasAuxInfo(layer_hash(cs)))
    weightvec =
        static_cast<weight_t*>(simulstate.getAuxInfo<void*>(layer_hash(cs)));

  if (!weightvec) {
    FCS_MSG_ERROR("Weights not stored in simulstate for calosample=" << cs);
    return FCSFatal;
  }

  TH2* voxel_template = nullptr;
  if (simulstate.hasAuxInfo(layer_hash_geo(cs)))
    voxel_template =
        static_cast<TH2*>(simulstate.getAuxInfo<void*>(layer_hash_geo(cs)));

  if (!voxel_template) {
    FCS_MSG_ERROR(
        "Voxel geometry not stored in simulstate for calosample=" << cs);
    return FCSFatal;
  }

  double deta_hit_minus_extrapol = hit.eta() - center_eta;
  const double dphi_hit_minus_extrapol =
      TVector2::Phi_mpi_pi(hit.phi() - center_phi);

  if (center_eta < 0.)
    deta_hit_minus_extrapol = -deta_hit_minus_extrapol;

  const double phi_dist2r = 1.0;
  const double dist000 = FastCaloSim::dist000(center_r, center_z);
  const double eta_jakobi = FastCaloSim::eta_jakobi(hit.eta());

  const double deta_hit_minus_extrapol_mm =
      deta_hit_minus_extrapol * eta_jakobi * dist000;
  const double dphi_hit_minus_extrapol_mm =
      dphi_hit_minus_extrapol * center_r * phi_dist2r;

  const double alpha_mm =
      std::atan2(dphi_hit_minus_extrapol_mm, deta_hit_minus_extrapol_mm);
  const double radius_mm =
      std::sqrt(dphi_hit_minus_extrapol_mm * dphi_hit_minus_extrapol_mm
                + deta_hit_minus_extrapol_mm * deta_hit_minus_extrapol_mm);

  int ix;
  int iy = voxel_template->GetYaxis()->FindBin(radius_mm) - 1;

  // Treat core as one bin for simulation
  if (iy == 0) {
    ix = 0;
  } else {
    if (alpha_mm < 0)
      ix = voxel_template->GetXaxis()->FindBin(2 * std::numbers::pi
                                               - std::fabs(alpha_mm))
          - 1;
    else
      ix = voxel_template->GetXaxis()->FindBin(alpha_mm) - 1;
  }

  double weight = 1;
  if (ix >= 0 && static_cast<std::size_t>(ix) < weightvec->size()) {
    const auto& weights_x = (*weightvec)[static_cast<std::size_t>(ix)];
    if (iy >= 0 && static_cast<std::size_t>(iy) < weights_x.size()) {
      weight = weights_x[static_cast<std::size_t>(iy)];
    }
  }

  double new_E = hit.E() * weight;
  hit.set_E(new_E);

  FCS_MSG_DEBUG("HIT: E=" << hit.E() << ", alpha = " << alpha_mm
                          << ", r = " << center_r << ", ix = " << ix
                          << ", iy = " << iy << ", weight = " << weight);

  return FCSSuccess;
}
