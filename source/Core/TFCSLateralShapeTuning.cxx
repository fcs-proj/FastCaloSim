// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#include <format>

#include "FastCaloSim/Core/TFCSLateralShapeTuning.h"

#include "FastCaloSim/Core/TFCSSimulationState.h"
#include "FastCaloSim/Definitions/ParticleData.h"
#include "FastCaloSim/Definitions/Phi.h"
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

auto TFCSLateralShapeTuning::initFromModelFile(
    const std::string& pathToModelParameters,
    int intMinEta,
    int intMaxEta,
    Model model) -> FCSReturnCode
{
  // get current calo layer
  std::size_t layer = TFCSLateralShapeParametrization::calosample();

  // load file containing model parameters
  std::unique_ptr<TFile> modelParametersFile = std::unique_ptr<TFile>(
      TFile::Open(pathToModelParameters.c_str(), "READ"));
  modelParametersFile->cd();

  // set parameter model names depending on the model
  const std::vector<std::string> parameterModelNames = model == Model::kSeries
      ? std::vector<std::string> {"a0", "a1", "a2", "a3"}
      : std::vector<std::string> {"eta_s", "phi_s"};

  for (auto& parameterName : parameterModelNames) {
    // interpolate the model parameters
    auto* linModelInterpol = new TFCSEnergyInterpolationPiecewiseLinear("", "");

    TGraph const* modelParameterGraph = static_cast<TGraph*>(
        gDirectory->Get(std::format("photons/layer{}/m{}_m{}_{}_{}/{}",
                                    layer,
                                    intMaxEta,
                                    intMinEta,
                                    intMinEta,
                                    intMaxEta,
                                    parameterName)
                            .c_str()));
    if (modelParameterGraph != nullptr) {
      // initialize the model interpolation and save in interpolation map
      linModelInterpol->InitFromArrayInEkin(modelParameterGraph->GetN(),
                                            modelParameterGraph->GetX(),
                                            modelParameterGraph->GetY());
      m_parameterInterpol.insert(
          std::make_pair(parameterName, linModelInterpol));
    } else {
      FCS_MSG_DEBUG(
          "[TFCSLateralShapeTuning] Could not find model parameter "
          "graph for layer="
          << layer << " minEta=" << intMinEta << " maxEta=" << intMaxEta);
      return FCSSuccess;
    }
  }
  modelParametersFile->Close();

  return FCSSuccess;
}

auto TFCSLateralShapeTuning::initFromMap(const interpolationMap& parameters)
    -> FCSReturnCode
{
  FCS_MSG_DEBUG(
      "[TFCSLateralShapeTuning] Initializing data tuning model from "
      "interpolation map.");
  m_parameterInterpol = parameters;
  return FCSSuccess;
}

auto TFCSLateralShapeTuning::simulate_hit(
    Hit& hit,
    TFCSSimulationState& /*simulstate*/,
    const TFCSTruthState* truth,
    const TFCSExtrapolationState* /*extrapol*/) -> FCSReturnCode
{
  // do not do anything if the parameter interpolation map is empty
  // this means we are in an pseudorapidity region, where no tuning to data is
  // available
  if (m_parameterInterpol.empty()) {
    return FCSSuccess;
  }

  // set maximum scaling factor
  const double maxScaling = 500;

  // retrieve particle data
  const int pdgId = truth->pdgid();
  const double charge = ParticleData::charge(pdgId);

  // retrieve hit information
  const double centerEta = hit.center_eta();
  const double centerPhi = hit.center_phi();
  const double centerZ = hit.center_z();
  // retrieve truth kinetic energy for interpolation
  const double Ekin = truth->Ekin();
  // retrieve calo sample
  std::size_t const layer = TFCSLateralShapeParametrization::calosample();

  FCS_MSG_DEBUG("[TFCSLateralShapeTuning] Initializing with pdgId="
                << pdgId << ", charge=" << charge << ", Ekin=" << Ekin
                << ", caloSample=" << layer);

  FCS_MSG_DEBUG("[TFCSLateralShapeTuning] Got hit position: "
                << " hit.eta=" << hit.eta() << ", hit.phi=" << hit.phi()
                << ", hit.r=" << hit.r());

  // compute deltaEta and deltaPhi
  const double deltaEta = hit.eta() - centerEta;
  const double deltaPhi = FastCaloSim::delta_phi(hit.phi(), centerPhi);

  // Apply the model whose parameters were loaded
  const auto etaScale = m_parameterInterpol.find("eta_s");
  if (etaScale != m_parameterInterpol.end()) {
    double etaScaleFactor = etaScale->second->evaluate(Ekin);
    double phiScaleFactor = m_parameterInterpol.at("phi_s")->evaluate(Ekin);

    // add a maximum scaling threshold to prevent unreasonable scalings
    etaScaleFactor =
        std::abs(etaScaleFactor) < maxScaling ? etaScaleFactor : maxScaling;
    phiScaleFactor =
        std::abs(phiScaleFactor) < maxScaling ? phiScaleFactor : maxScaling;

    FCS_MSG_DEBUG(
        "[TFCSLateralShapeTuning] Applying 2D eta_s - eta_phi "
        "scaling model with eta_s="
        << etaScaleFactor << " and phi_s=" << phiScaleFactor);

    double const deltaEtaCorr = etaScaleFactor * deltaEta;
    double const deltaPhiCorr = phiScaleFactor * deltaPhi;

    hit.setEtaPhiZE(centerEta + deltaEtaCorr,
                    FastCaloSim::wrap_phi(centerPhi + deltaPhiCorr),
                    centerZ,
                    hit.E());

  } else {
    double etaScaleFactor =
        getSeriesScalingFactor(m_parameterInterpol.at("a0")->evaluate(Ekin),
                               m_parameterInterpol.at("a1")->evaluate(Ekin),
                               m_parameterInterpol.at("a2")->evaluate(Ekin),
                               m_parameterInterpol.at("a3")->evaluate(Ekin),
                               std::abs(deltaEta));

    // add a maximum scaling threshold to prevent unreasonable scalings
    etaScaleFactor =
        std::abs(etaScaleFactor) < maxScaling ? etaScaleFactor : maxScaling;

    FCS_MSG_DEBUG(
        "[TFCSLateralShapeTuning] Applying eta_s series expansion "
        "model with eta_s="
        << etaScaleFactor);

    double const deltaEtaCorr = etaScaleFactor * deltaEta;

    hit.setEtaPhiZE(centerEta + deltaEtaCorr,
                    FastCaloSim::wrap_phi(centerPhi + deltaPhi),
                    centerZ,
                    hit.E());
  }

  return FCSSuccess;
}

auto TFCSLateralShapeTuning::getSeriesScalingFactor(
    double a0, double a1, double a2, double a3, double distToShowerCenter)
    -> double
{
  double const meanDistToShowerCentre = 0.0039;
  double const scaleFactor = 1 + a0
      + a1 * distToShowerCenter / meanDistToShowerCentre
      + a2 * std::pow(distToShowerCenter / meanDistToShowerCentre, 2)
      + a3 * std::pow(distToShowerCenter / meanDistToShowerCentre, 3);

  return scaleFactor;
}
