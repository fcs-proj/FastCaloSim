/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#include "FastCaloSim/Core/TFCSLateralShapeTuning.h"

#include "FastCaloSim/Core/TFCSSimulationState.h"
#include "FastCaloSim/Definitions/ParticleData.h"
#include "TFile.h"
#include "TGraph.h"

//=============================================
//======= TFCSLateralShapeTuning ==============
//=============================================

TFCSLateralShapeTuning::TFCSLateralShapeTuning(const char* name,
                                               const char* title)
    : TFCSLateralShapeParametrizationHitBase(name, title)
{
}

TFCSLateralShapeTuning::~TFCSLateralShapeTuning()
{
  // clear parameter interpolation map and free memory
  for (auto& p : m_parameterInterpol) {
    delete p.second;
  }
  m_parameterInterpol.clear();
}

FCSReturnCode TFCSLateralShapeTuning::initFromModelFile(
    const std::string& pathToModelParameters, int intMinEta, int intMaxEta)
{
  // get current calo layer
  int layer = TFCSLateralShapeParametrization::calosample();

  // load file containing model parameters
  std::unique_ptr<TFile> modelParametersFile = std::unique_ptr<TFile>(
      TFile::Open(pathToModelParameters.c_str(), "READ"));
  modelParametersFile->cd();

  // set parameter model names depending on layer
  std::vector<std::string> parameterModelNames;
  if (layer == 1 || layer == 5)
    parameterModelNames = {"a0", "a1", "a2", "a3"};
  else if (layer == 2 || layer == 6)
    parameterModelNames = {"eta_s", "phi_s"};

  for (auto& parameterName : parameterModelNames) {
    // interpolate the model parameters
    TFCSEnergyInterpolationPiecewiseLinear* linModelInterpol =
        new TFCSEnergyInterpolationPiecewiseLinear("", "");

    TGraph* modelParameterGraph =
        (TGraph*)gDirectory->Get(Form("photons/layer%d/m%d_m%d_%d_%d/%s",
                                      layer,
                                      intMaxEta,
                                      intMinEta,
                                      intMinEta,
                                      intMaxEta,
                                      parameterName.c_str()));
    if (modelParameterGraph) {
      // initialize the model interpolation and save in interpolation map
      linModelInterpol->InitFromArrayInEkin(modelParameterGraph->GetN(),
                                            modelParameterGraph->GetX(),
                                            modelParameterGraph->GetY());
      m_parameterInterpol.insert(
          std::make_pair(parameterName, linModelInterpol));
    } else {
      ATH_MSG_DEBUG(
          "[TFCSLateralShapeTuning] Could not find model parameter "
          "graph for layer="
          << layer << " minEta=" << intMinEta << " maxEta=" << intMaxEta);
      return FCSSuccess;
    }
  }
  modelParametersFile->Close();

  return FCSSuccess;
}

FCSReturnCode TFCSLateralShapeTuning::initFromMap(
    const interpolationMap& interpolationMap)
{
  ATH_MSG_DEBUG(
      "[TFCSLateralShapeTuning] Initializing data tuning model from "
      "interpolation map.");
  m_parameterInterpol = interpolationMap;
  return FCSSuccess;
}

FCSReturnCode TFCSLateralShapeTuning::simulate_hit(
    Hit& hit,
    TFCSSimulationState&,
    const TFCSTruthState* truth,
    const TFCSExtrapolationState*)
{
  // do not do anything if the parameter interpolation map is empty
  // this means we are in an pseudorapidity region, where no tuning to data is
  // available
  if (m_parameterInterpol.empty())
    return FCSSuccess;

  // set maximum scaling factor
  float maxScaling = 500;

  // retrieve particle data
  const int pdgId = truth->pdgid();
  const double charge = ParticleData::charge(pdgId);

  // retrieve hit information
  const double centerEta = hit.center_eta();
  const double centerPhi = hit.center_phi();
  const double centerZ = hit.center_z();
  // retrieve truth kinetic energy for interpolation
  const float Ekin = truth->Ekin();
  // retrieve calo sample
  int layer = TFCSLateralShapeParametrization::calosample();

  ATH_MSG_DEBUG("[TFCSLateralShapeTuning] Initializing with pdgId="
                << pdgId << ", charge=" << charge << ", Ekin=" << Ekin
                << ", caloSample=" << layer);

  ATH_MSG_DEBUG("[TFCSLateralShapeTuning] Got hit position: "
                << " hit.eta=" << hit.eta() << ", hit.phi=" << hit.phi()
                << ", hit.r=" << hit.r());

  // compute deltaEta and deltaPhi
  const double deltaEta = hit.eta() - centerEta;
  const double deltaPhi = hit.phi() - centerPhi;

  if (layer == 2 || layer == 6) {
    double etaScaleFactor = m_parameterInterpol["eta_s"]->evaluate(Ekin);
    double phiScaleFactor = m_parameterInterpol["phi_s"]->evaluate(Ekin);

    // add a maximum scaling threshold to prevent unreasonable scalings
    etaScaleFactor =
        std::abs(etaScaleFactor) < maxScaling ? etaScaleFactor : maxScaling;
    phiScaleFactor =
        std::abs(phiScaleFactor) < maxScaling ? phiScaleFactor : maxScaling;

    ATH_MSG_DEBUG(
        "[TFCSLateralShapeTuning] Applying 2D eta_s - eta_phi "
        "scaling model with eta_s="
        << etaScaleFactor << " and phi_s=" << phiScaleFactor);

    double deltaEtaCorr = etaScaleFactor * deltaEta;
    double deltaPhiCorr = phiScaleFactor * deltaPhi;

    hit.setEtaPhiZE(
        centerEta + deltaEtaCorr, centerPhi + deltaPhiCorr, centerZ, hit.E());

  }

  else if (layer == 1 || layer == 5)
  {
    double etaScaleFactor =
        getSeriesScalingFactor(m_parameterInterpol["a0"]->evaluate(Ekin),
                               m_parameterInterpol["a1"]->evaluate(Ekin),
                               m_parameterInterpol["a2"]->evaluate(Ekin),
                               m_parameterInterpol["a3"]->evaluate(Ekin),
                               std::abs(deltaEta));

    // add a maximum scaling threshold to prevent unreasonable scalings
    etaScaleFactor =
        std::abs(etaScaleFactor) < maxScaling ? etaScaleFactor : maxScaling;

    ATH_MSG_DEBUG(
        "[TFCSLateralShapeTuning] Applying eta_s series expansion "
        "model with eta_s="
        << etaScaleFactor);

    double deltaEtaCorr = etaScaleFactor * deltaEta;

    hit.setEtaPhiZE(
        centerEta + deltaEtaCorr, centerPhi + deltaPhi, centerZ, hit.E());
  }

  return FCSSuccess;
}

double TFCSLateralShapeTuning::getSeriesScalingFactor(
    double a0, double a1, double a2, double a3, double distToShowerCenter)
{
  double meanDistToShowerCentre = 0.0039;
  double scaleFactor = 1 + a0 + a1 * distToShowerCenter / meanDistToShowerCentre
      + a2 * std::pow(distToShowerCenter / meanDistToShowerCentre, 2)
      + a3 * std::pow(distToShowerCenter / meanDistToShowerCentre, 3);

  return scaleFactor;
}
