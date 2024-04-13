/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#include <string>

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGauss.h"

#include "FastCaloSim/TFCSVoxelHistoLateralCovarianceFluctuations.h"
#include "FastCaloSim/ICaloGeometry.h"
#include "FastCaloSim/TFCSSimulationState.h"
#include "FastCaloSim/TFCS1DFunction.h"
#include "FastCaloSim/TFCS1DFunctionTemplateHistogram.h"
#include "FastCaloSim/ICaloGeometry.h"

#include "TVector2.h"
#include "TMath.h"
#include "TFile.h"
#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"
#include "TH2.h"

//=============================================
//======= TFCSVoxelHistoLateralCovarianceFluctuations =========
//=============================================

const std::uint32_t TFCSVoxelHistoLateralCovarianceFluctuations::s_layer_hash
    [CaloCell_ID_FCS::MaxSample] = {
        "PreSamplerB"_FCShash, "EMB1"_FCShash,        "EMB2"_FCShash,
        "EMB3"_FCShash,        "PreSamplerE"_FCShash, "EME1"_FCShash,
        "EME2"_FCShash,        "EME3"_FCShash,        "HEC0"_FCShash,
        "HEC1"_FCShash,        "HEC2"_FCShash,        "HEC3"_FCShash,
        "TileBar0"_FCShash,    "TileBar1"_FCShash,    "TileBar2"_FCShash,
        "TileGap1"_FCShash,    "TileGap2"_FCShash,    "TileGap3"_FCShash,
        "TileExt0"_FCShash,    "TileExt1"_FCShash,    "TileExt2"_FCShash,
        "FCAL0"_FCShash,       "FCAL1"_FCShash,       "FCAL2"_FCShash};

const std::uint32_t TFCSVoxelHistoLateralCovarianceFluctuations::
    s_layer_hash_geo[CaloCell_ID_FCS::MaxSample] = {
        "PreSamplerB_geo"_FCShash, "EMB1_geo"_FCShash,
        "EMB2_geo"_FCShash,        "EMB3_geo"_FCShash,
        "PreSamplerE_geo"_FCShash, "EME1_geo"_FCShash,
        "EME2_geo"_FCShash,        "EME3_geo"_FCShash,
        "HEC0_geo"_FCShash,        "HEC1_geo"_FCShash,
        "HEC2_geo"_FCShash,        "HEC3_geo"_FCShash,
        "TileBar0_geo"_FCShash,    "TileBar1_geo"_FCShash,
        "TileBar2_geo"_FCShash,    "TileGap1_geo"_FCShash,
        "TileGap2_geo"_FCShash,    "TileGap3_geo"_FCShash,
        "TileExt0_geo"_FCShash,    "TileExt1_geo"_FCShash,
        "TileExt2_geo"_FCShash,    "FCAL0_geo"_FCShash,
        "FCAL1_geo"_FCShash,       "FCAL2_geo"_FCShash};

TFCSVoxelHistoLateralCovarianceFluctuations::
    TFCSVoxelHistoLateralCovarianceFluctuations(const char *name,
                                                const char *title)
    : TFCSLateralShapeParametrizationHitBase(name, title) {}

TFCSVoxelHistoLateralCovarianceFluctuations::
    ~TFCSVoxelHistoLateralCovarianceFluctuations() {}

bool TFCSVoxelHistoLateralCovarianceFluctuations::initialize(
    TFile *inputfile, const std::string &folder) {
  // load m_eigenvariances and m_parMeans from input file
  // load histograms for each cell from input file
  // m_transform[i][j]=new TFCS1DFunctionInt16Int16Histogram(hist)
  // where hist is the 1D histogram from which a CDF was calculated

  inputfile->cd(folder.c_str());
  //  inputfile->ls();

  int cs = calosample();
  int bin = 1;

  TH2 *temp = dynamic_cast<TH2 *>(
      inputfile->Get(Form("voxel_template_cs%d_pca%d", cs, bin)));
  if (!temp) {
    ATH_MSG_ERROR("Template hist not found for cs " + std::to_string(cs));
    return false;
  }
  m_voxel_template.push_back(temp);
  m_nDim_x = m_voxel_template[0]->GetNbinsX();
  m_nDim_y = m_voxel_template[0]->GetNbinsY();

  while (inputfile->Get(Form("parMeans_cs%d_pca%d", cs, bin))) {
    TVectorD parMeans;
    TMatrixD EigenVectors;
    TVectorD EigenValues;
    std::vector<std::vector<TFCS1DFunction *>> transform;

    std::string label = Form("_cs%d_pca%d", cs, bin);

    TObject *obj = inputfile->Get(("parMeans" + label).c_str());
    if (!obj) {
      ATH_MSG_ERROR("parMeans" + label + " not found");
      return false;
    }
    if (msgLvl(MSG::DEBUG)) {
      ATH_MSG_DEBUG("parMeans");
      obj->Print();
    }
    parMeans.ResizeTo(m_nDim_x * m_nDim_y);
    parMeans = *dynamic_cast<TVectorT<double> *>(obj);
    m_parMeans.push_back(parMeans);

    obj = inputfile->Get(("covMatrix" + label).c_str());
    if (!obj) {
      ATH_MSG_ERROR("covMatrix" + label + " not found");
      return false;
    }
    TMatrixTSym<double> covMatrix = *dynamic_cast<TMatrixTSym<double> *>(obj);
    TMatrixDSymEigen eigenvariances = TMatrixDSymEigen(covMatrix);
    EigenVectors.ResizeTo(eigenvariances.GetEigenVectors());
    EigenVectors = eigenvariances.GetEigenVectors();
    m_EigenVectors.push_back(EigenVectors);

    EigenValues.ResizeTo(eigenvariances.GetEigenValues());
    EigenValues = eigenvariances.GetEigenValues();
    m_EigenValues.push_back(EigenValues);

    if (msgLvl(MSG::DEBUG)) {
      ATH_MSG_DEBUG("eigenvariances");
      eigenvariances.GetEigenValues().Print();
      eigenvariances.GetEigenVectors().Print();
    }
    if (msgLvl(MSG::DEBUG)) {
      ATH_MSG_DEBUG("m_EigenValues");
      EigenValues.Print();
      ATH_MSG_DEBUG("m_EigenVectors");
      EigenVectors.Print();
    }

    transform.resize(m_nDim_x);
    for (int x = 0; x < m_nDim_x; ++x) {
      transform[x].resize(m_nDim_y);
      for (int y = 0; y < m_nDim_y; ++y) {
        // Naming of hists is hist_xy where x, y are indices of cells
        std::string histname = Form("hist_%d_%d%s", x, y, label.c_str());
        TH1 *hist = (TH1 *)inputfile->Get(histname.c_str());
        if (!hist) {
          ATH_MSG_ERROR("Histogram " << histname << " not found");
          return false;
        }
        TFCS1DFunctionInt32Int32Histogram *func =
            new TFCS1DFunctionInt32Int32Histogram();
        func->Initialize(hist, false);
        transform[x][y] = func;
      }
    }
    m_transform.push_back(transform);
    ++bin;
  }

  return true;
}

// Implementation of multivariate Gaussian with ROOT classes (from ROOT forum,
// but validated)
void TFCSVoxelHistoLateralCovarianceFluctuations::MultiGaus(
    TFCSSimulationState &simulstate, TVectorD &genPars) const {
  int Ebin = simulstate.Ebin();

  int nPars = m_parMeans[Ebin - 1].GetNrows();
  genPars.ResizeTo(nPars);

  TVectorD rotParMeans = m_EigenVectors[Ebin - 1] * m_parMeans[Ebin - 1];
  for (int iPar = 0; iPar < nPars; iPar++) {
    double variance = m_EigenValues[Ebin - 1][iPar];
    // check for positive-definiteness of covMatrix
    if (variance < 0) {
      ATH_MSG_ERROR("Got a negative eigenvariance (" << variance
                                                     << ") on iPar = " << iPar);
      variance = 0;
    }
    genPars[iPar] = CLHEP::RandGauss::shoot(simulstate.randomEngine(),
                                            rotParMeans[iPar], sqrt(variance));
    ATH_MSG_DEBUG("genPars[" << iPar << "]=" << genPars[iPar]
                             << " rotParMeans[iPar]=" << rotParMeans[iPar]
                             << " sqrt(variance)=" << sqrt(variance));
  }
  genPars = m_EigenVectors[Ebin - 1] * genPars;
}

FCSReturnCode TFCSVoxelHistoLateralCovarianceFluctuations::simulate(
    TFCSSimulationState &simulstate, const TFCSTruthState * /*truth*/,
    const TFCSExtrapolationState * /*extrapol*/) const {
  if (!simulstate.randomEngine()) {
    return FCSFatal;
  }

  if (msgLvl(MSG::DEBUG)) {
    ATH_MSG_DEBUG("simulstate before clearing AuxInfo");
    simulstate.Print();
  }

  weight_t *weightvec;
  TH2 *voxel_temp;

  int Ebin = simulstate.Ebin();

  for (int ilayer = 0; ilayer < CaloCell_ID_FCS::MaxSample; ++ilayer) {
    if (simulstate.hasAuxInfo(s_layer_hash[ilayer])) {
      weightvec = static_cast<weight_t *>(
          simulstate.getAuxInfo<void *>(s_layer_hash[ilayer]));
      if (weightvec) {
        delete weightvec;
        simulstate.setAuxInfo<void *>(s_layer_hash[ilayer], nullptr);
      }
    }
    if (simulstate.hasAuxInfo(s_layer_hash_geo[ilayer])) {
      voxel_temp = static_cast<TH2 *>(
          simulstate.getAuxInfo<void *>(s_layer_hash_geo[ilayer]));
      if (voxel_temp) {
        delete voxel_temp;
        simulstate.setAuxInfo<void *>(s_layer_hash_geo[ilayer], nullptr);
      }
    }
  }

  if (msgLvl(MSG::DEBUG)) {
    ATH_MSG_DEBUG("simulstate after clearing AuxInfo");
    simulstate.Print();
  }
  // if(m_parMeans[Ebin-1].GetNrows()<1) {
  if (m_parMeans.empty() || Ebin <= 0) {
    // not initialized, do nothing in case this instance is just used to clean
    // up memory
    return FCSSuccess;
  }

  // TODO: the following code should be executed for all relevant calo layers,
  // possibly simulating correlated fluctuations between layers depending on the
  // PCA bin
  weightvec = new weight_t(m_nDim_x);
  for (int x = 0; x < m_nDim_x; ++x)
    (*weightvec)[x].resize(m_nDim_y);

  TVectorD genPars(m_parMeans[Ebin - 1].GetNrows());

  // Fill genPars with multivarate random Gauss
  MultiGaus(simulstate, genPars);

  // Loop through voxels
  int count = 0;
  for (int x = 0; x < m_nDim_x; ++x) {
    for (int y = 0; y < m_nDim_y; ++y) {
      // Get value of CDF corresponding to generated Gauss random
      double cdf_val =
          (TMath::Erf(genPars[count] * 2.0 / TMath::Pi()) + 1) / 2.0;

      // Map that cdf value to input bin and value
      double orig_val = m_transform[Ebin - 1][x][y]->rnd_to_fct(cdf_val);

      // Fill hist, keep linear count
      (*weightvec)[x][y] = orig_val;
      ++count;

      ATH_MSG_DEBUG("CELL x=" << x << " y=" << y << " : cdf_val=" << cdf_val
                              << " orig_val=" << orig_val);
    }
  }

  voxel_temp = static_cast<TH2 *>(m_voxel_template[0]->Clone());

  // For now simulating only the layer calosample()
  simulstate.setAuxInfo<void *>(s_layer_hash[calosample()], weightvec);
  simulstate.setAuxInfo<void *>(s_layer_hash_geo[calosample()], voxel_temp);

  if (msgLvl(MSG::DEBUG)) {
    ATH_MSG_DEBUG("simulstate after storing weight " << weightvec
                                                     << " in AuxInfo");
    simulstate.Print();
  }

  return FCSSuccess;
}

FCSReturnCode TFCSVoxelHistoLateralCovarianceFluctuations::simulate_hit(
    Hit &hit, TFCSSimulationState &simulstate, const TFCSTruthState * /*truth*/,
    const TFCSExtrapolationState * /*extrapol*/) {
  const double center_eta = hit.center_eta();
  const double center_phi = hit.center_phi();
  const double center_r = hit.center_r();
  const double center_z = hit.center_z();

  int cs = calosample();
  weight_t *weightvec = nullptr;
  if (simulstate.hasAuxInfo(s_layer_hash[cs]))
    weightvec = static_cast<weight_t *>(
        simulstate.getAuxInfo<void *>(s_layer_hash[cs]));

  if (!weightvec) {
    ATH_MSG_ERROR("Weights not stored in simulstate for calosample=" << cs);
    return FCSFatal;
  }

  TH2 *voxel_template = nullptr;
  if (simulstate.hasAuxInfo(s_layer_hash_geo[cs]))
    voxel_template =
        static_cast<TH2 *>(simulstate.getAuxInfo<void *>(s_layer_hash_geo[cs]));
  if (!voxel_template) {
    ATH_MSG_ERROR(
        "Voxel geometry not stored in simulstate for calosample=" << cs);
    return FCSFatal;
  }

  /*  int nbinsx = 8;
    float lowx = 0;
    float highx = 2 * TMath::Pi();
    int nbinsy = 9.;
    float lowy = 0.;
    float highy = 165.;
    float coreSize = 5.;

    std::vector<float> bins_x = {};
    for(int bin = 0; bin <= nbinsx; ++bin){
        bins_x.push_back(lowx+((highx-lowx)/nbinsx)*bin);
    }

    std::vector<float> bins_y = {lowy};
    for(int bin = 0; bin <= nbinsy; ++bin){
        float step = (highy-lowy)/nbinsy;
        if(coreSize > 0){
            if(bin ==0){lowy =lowy+coreSize;}
            bins_y.push_back(lowy+bin*step);
        }
        else{
            if(bin > 0){
                bins_y.push_back(lowy+bin*step);
            }
        }
    }

    TH2F* voxel_template = new TH2F("", "", bins_x.size()-1, &bins_x[0],
    bins_y.size()-1, &bins_y[0]); */

  float deta_hit_minus_extrapol = hit.eta() - center_eta;
  float dphi_hit_minus_extrapol = TVector2::Phi_mpi_pi(hit.phi() - center_phi);

  //  if(charge<0.)dphi_hit_minus_extrapol=-dphi_hit_minus_extrapol;
  if (center_eta < 0.)
    deta_hit_minus_extrapol = -deta_hit_minus_extrapol;

  float phi_dist2r = 1.0;
  float dist000 = TMath::Sqrt(center_r * center_r + center_z * center_z);
  float eta_jakobi = TMath::Abs(2.0 * TMath::Exp(-hit.eta()) /
                                (1.0 + TMath::Exp(-2 * hit.eta())));

  float deta_hit_minus_extrapol_mm =
      deta_hit_minus_extrapol * eta_jakobi * dist000;
  float dphi_hit_minus_extrapol_mm =
      dphi_hit_minus_extrapol * center_r * phi_dist2r;

  float alpha_mm =
      TMath::ATan2(dphi_hit_minus_extrapol_mm, deta_hit_minus_extrapol_mm);
  float radius_mm =
      TMath::Sqrt(dphi_hit_minus_extrapol_mm * dphi_hit_minus_extrapol_mm +
                  deta_hit_minus_extrapol_mm * deta_hit_minus_extrapol_mm);

  int ix;
  int iy = voxel_template->GetYaxis()->FindBin(radius_mm) - 1;

  // Treat core as one bin for simulation
  if (iy == 0) {
    ix = 0;
  } else {
    if (alpha_mm < 0)
      ix = voxel_template->GetXaxis()->FindBin(2 * TMath::Pi() -
                                               fabs(alpha_mm)) -
           1;
    else
      ix = voxel_template->GetXaxis()->FindBin(alpha_mm) - 1;
  }

  const int sizex = (*weightvec).size();

  float weight = 1;
  if (ix >= 0 && ix < sizex) {
    const int sizey = (*weightvec)[ix].size();
    if (iy >= 0 && iy < sizey) {
      weight = (*weightvec)[ix][iy];
    }
  }

  hit.E() *= weight;

  ATH_MSG_DEBUG("HIT: E=" << hit.E() << ", alpha = " << alpha_mm
                          << ", r = " << center_r << ", ix = " << ix
                          << ", iy = " << iy << ", weight = " << weight);

  // delete voxel_template;
  return FCSSuccess;
}
