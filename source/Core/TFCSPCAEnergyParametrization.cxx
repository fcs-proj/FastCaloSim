// Copyright (c) 2024 CERN for the benefit of the FastCaloSim project

#include "FastCaloSim/Core/TFCSPCAEnergyParametrization.h"

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGaussZiggurat.h"
#include "FastCaloSim/Core/TFCSExtrapolationState.h"
#include "FastCaloSim/Core/TFCSSimulationState.h"
#include "FastCaloSim/Geometry/CaloGeo.h"
#include "TClass.h"
#include "TFile.h"
#include "TH1.h"
#include "TKey.h"
#include "TMath.h"
#include "TMatrixD.h"
#include "TMatrixDSymEigen.h"

//=============================================
//======= TFCSPCAEnergyParametrization =========
//=============================================

TFCSPCAEnergyParametrization::TFCSPCAEnergyParametrization(const char* name,
                                                           const char* title,
                                                           CaloGeo* geo)
    : TFCSEnergyParametrization(name, title)
    , m_geo(geo)
{
  m_numberpcabins = 1;
  do_rescale = 1;
}

bool TFCSPCAEnergyParametrization::is_match_Ekin_bin(int Ekin_bin) const
{
  if (Ekin_bin >= 1 && Ekin_bin <= n_bins())
    return true;
  return false;
}

bool TFCSPCAEnergyParametrization::is_match_calosample(int calosample) const
{
  for (unsigned int i = 0; i < m_RelevantLayers.size(); i++) {
    if (m_RelevantLayers[i] == calosample)
      return true;
  }
  return false;
}

void TFCSPCAEnergyParametrization::Print(Option_t* option) const
{
  TString opt(option);
  bool shortprint = opt.Index("short") >= 0;
  bool longprint =
      msgLvl(FCS_MSG::DEBUG) || (msgLvl(FCS_MSG::INFO) && !shortprint);
  TString optprint = opt;
  optprint.ReplaceAll("short", "");
  TFCSEnergyParametrization::Print(option);

  if (longprint) {
    FCS_MSG(INFO) << optprint << "  #bins=" << m_numberpcabins
                  << ", Enorm=" << m_total_energy_normalization << ", layers=";
    for (unsigned int i = 0; i < m_RelevantLayers.size(); i++) {
      if (i > 0)
        FCS_MSG(INFO) << ", ";
      FCS_MSG(INFO) << m_RelevantLayers[i];
    }
    FCS_MSG(INFO) << END_FCS_MSG(INFO);
  }
}

float interpolate_get_y(TH1* hist, float x)
{
  float m = 0;
  float n = 0;
  float x1 = 1;
  float x2 = 0;
  float y1 = 0;
  float y2 = 0;
  if (x <= hist->GetBinCenter(1))
    return hist->GetBinContent(1);
  if (x >= hist->GetBinCenter(hist->GetNbinsX()))
    return hist->GetBinContent(hist->GetNbinsX());

  int bin = hist->FindBin(x);
  // is x above the bin center -> interpolate with next bin
  if (x > hist->GetBinCenter(bin)) {
    x1 = hist->GetBinCenter(bin);
    y1 = hist->GetBinContent(bin);
    x2 = hist->GetBinCenter(bin + 1);
    y2 = hist->GetBinContent(bin + 1);
  }

  // is x below bin center -> interpolate with previous bin
  if (x <= hist->GetBinCenter(bin)) {
    x1 = hist->GetBinCenter(bin - 1);
    y1 = hist->GetBinContent(bin - 1);
    x2 = hist->GetBinCenter(bin);
    y2 = hist->GetBinContent(bin);
  }

  m = (y1 - y2) / (x1 - x2);
  n = y2 - m * x2;
  return m * x + n;
}

void TFCSPCAEnergyParametrization::set_totalE_probability_ratio(int Ekin_bin,
                                                                TH1* hist)
{
  if (Ekin_bin < 1)
    return;
  if (Ekin_bin - 1 >= (int)m_totalE_probability_ratio.size())
    return;
  m_totalE_probability_ratio[Ekin_bin - 1] = hist;
}

TH1* TFCSPCAEnergyParametrization::get_totalE_probability_ratio(
    int Ekin_bin) const
{
  if (Ekin_bin < 1)
    return nullptr;
  if (Ekin_bin - 1 >= (int)m_totalE_probability_ratio.size())
    return nullptr;
  return m_totalE_probability_ratio[Ekin_bin - 1];
}

FCSReturnCode TFCSPCAEnergyParametrization::simulate(
    TFCSSimulationState& simulstate,
    const TFCSTruthState* /*truth*/,
    const TFCSExtrapolationState* /*extrapol*/) const
{
  if (!simulstate.randomEngine()) {
    return FCSFatal;
  }

  int pcabin = simulstate.Ebin();

  if (pcabin == 0) {
    simulstate.set_E(0);
    for (int s = 0; s < m_geo->n_layers(); s++) {
      simulstate.set_E(s, 0.0);
      simulstate.set_Efrac(s, 0.0);
    }
  } else {
    TMatrixD* EV = m_EV[pcabin - 1];
    TVectorD* MeanValues = m_MeanValues[pcabin - 1];
    TVectorD* SigmaValues = m_SigmaValues[pcabin - 1];
    TVectorD* Gauss_means = m_Gauss_means[pcabin - 1];
    TVectorD* Gauss_rms = m_Gauss_rms[pcabin - 1];
    std::vector<TFCS1DFunction*> cumulative = m_cumulative[pcabin - 1];

    const std::vector<int>& layerNr = m_RelevantLayers;

    double* vals_gauss_means = (double*)Gauss_means->GetMatrixArray();
    double* vals_gauss_rms = Gauss_rms->GetMatrixArray();

    double* output_data = new double[layerNr.size() + 1];
    double* input_data = new double[layerNr.size() + 1];

    for (unsigned int l = 0; l <= layerNr.size(); l++) {
      double mean = vals_gauss_means[l];
      double rms = vals_gauss_rms[l];
      double gauszz =
          CLHEP::RandGaussZiggurat::shoot(simulstate.randomEngine(), mean, rms);
      input_data[l] = gauszz;
    }

    P2X(SigmaValues,
        MeanValues,
        EV,
        layerNr.size() + 1,
        input_data,
        output_data,
        layerNr.size() + 1);

    double* simdata = new double[layerNr.size() + 1];
    double sum_fraction = 0.0;
    for (unsigned int l = 0; l <= layerNr.size(); l++) {
      double simdata_uniform =
          (TMath::Erf(output_data[l] / 1.414213562) + 1) / 2.f;

      simdata[l] = cumulative[l]->rnd_to_fct(simdata_uniform);

      if (l != layerNr.size())  // sum up the fractions, but not the totalE
        sum_fraction += simdata[l];
    }

    double scalefactor = 1.0 / sum_fraction;
    if (!do_rescale)
      scalefactor = 1.0;

    // Using signed int for better loop optimization and vectorization where
    // possible
    int layerNrsize = static_cast<int>(layerNr.size());
    for (int l = 0; l < layerNrsize; l++) {
      simdata[l] *= scalefactor;
    }

    // Apply hit-and-miss reweighting of total energy response
    // This needs to run before simulstate is modified, otherwise a clean retry
    // is not possible
    TH1* h_totalE_ratio = get_totalE_probability_ratio(pcabin);
    if (h_totalE_ratio) {
      float pass_probability = 0;
      float Etot = simdata[layerNr.size()];
      if (Etot > h_totalE_ratio->GetXaxis()->GetXmin()
          && Etot < h_totalE_ratio->GetXaxis()->GetXmax())
      {
        pass_probability = interpolate_get_y(h_totalE_ratio, Etot);
      }
      float random = CLHEP::RandFlat::shoot(simulstate.randomEngine());
      if (random > pass_probability) {
        delete[] output_data;
        delete[] input_data;
        delete[] simdata;

        return (FCSReturnCode)FCSRetryPCA;
      }
    }

    double total_energy =
        simdata[layerNr.size()] * simulstate.E() / m_total_energy_normalization;
    simulstate.set_E(total_energy);
    FCS_MSG_DEBUG("set E to total_energy=" << total_energy);

    for (int s = 0; s < m_geo->n_layers(); s++) {
      double energyfrac = 0.0;
      for (unsigned int l = 0; l < layerNr.size(); l++) {
        if (layerNr[l] == s)
          energyfrac = simdata[l];
      }
      simulstate.set_Efrac(s, energyfrac);
      simulstate.set_E(s, energyfrac * total_energy);
      simulstate.set_SF(scalefactor);
    }

    delete[] output_data;
    delete[] input_data;
    delete[] simdata;
  }

  return FCSSuccess;
}

void TFCSPCAEnergyParametrization::P2X(TVectorD* SigmaValues,
                                       TVectorD* MeanValues,
                                       TMatrixD* EV,
                                       int gNVariables,
                                       const double* p,
                                       double* x,
                                       int nTest)
{
  const double* gSigmaValues = SigmaValues->GetMatrixArray();
  const double* gMeanValues = MeanValues->GetMatrixArray();
  const double* gEigenVectors = EV->GetMatrixArray();

  for (int i = 0; i < gNVariables; i++) {
    x[i] = gMeanValues[i];
    for (int j = 0; j < nTest; j++) {
      x[i] +=
          p[j] * gSigmaValues[i] * (double)(gEigenVectors[i * gNVariables + j]);
    }
  }
}

bool TFCSPCAEnergyParametrization::loadInputs(TFile* file)
{
  return loadInputs(file, "");
}

bool TFCSPCAEnergyParametrization::loadInputs(TFile* file,
                                              const std::string& folder)
{
  bool load_ok = 1;

  int trynext = 1;
  TString x;
  if (folder.empty())
    x = "bin";
  else
    x = folder + "/bin";
  while (trynext) {
    IntArray* test = (IntArray*)file->Get(
        x + Form("%i/pca/RelevantLayers", m_numberpcabins));
    if (test) {
      m_numberpcabins++;
      delete test;
    } else
      trynext = 0;
  }
  m_numberpcabins -= 1;

  file->cd(x + "1/pca");
  IntArray* RelevantLayers = (IntArray*)gDirectory->Get("RelevantLayers");
  if (RelevantLayers == nullptr) {
    FCS_MSG_ERROR(
        "TFCSPCAEnergyParametrization::m_RelevantLayers in first "
        "pcabin is null!");
    load_ok = false;
  }

  if (!load_ok)
    return false;

  m_RelevantLayers.reserve(RelevantLayers->GetSize());
  for (int i = 0; i < RelevantLayers->GetSize(); i++)
    m_RelevantLayers.push_back(RelevantLayers->GetAt(i));

  for (int bin = 1; bin <= m_numberpcabins; bin++) {
    file->cd(x + Form("%i/pca", bin));

    TMatrixDSym* symCov = (TMatrixDSym*)gDirectory->Get("symCov");
    TVectorD* MeanValues = (TVectorD*)gDirectory->Get("MeanValues");
    TVectorD* SigmaValues = (TVectorD*)gDirectory->Get("SigmaValues");
    TVectorD* Gauss_means = (TVectorD*)gDirectory->Get("Gauss_means");
    TVectorD* Gauss_rms = (TVectorD*)gDirectory->Get("Gauss_rms");

    if (symCov == nullptr) {
      FCS_MSG_WARNING("TFCSPCAEnergyParametrization::symCov in pcabin "
                      << bin << " is null!");
      load_ok = false;
    }
    if (MeanValues == nullptr) {
      FCS_MSG_WARNING("TFCSPCAEnergyParametrization::MeanValues in pcabin "
                      << bin << " is null!");
      load_ok = false;
    }
    if (SigmaValues == nullptr) {
      FCS_MSG_WARNING("TFCSPCAEnergyParametrization::SigmaValues in pcabin "
                      << bin << " is null!");
      load_ok = false;
    }
    if (Gauss_means == nullptr) {
      FCS_MSG_WARNING("TFCSPCAEnergyParametrization::Gauss_means in pcabin "
                      << bin << " is null!");
      load_ok = false;
    }
    if (Gauss_rms == nullptr) {
      FCS_MSG_WARNING("TFCSPCAEnergyParametrization::Gause_rms in pcabin "
                      << bin << " is null!");
      load_ok = false;
    }

    if (!load_ok)
      return false;

    TMatrixDSymEigen cov_eigen(*symCov);
    TMatrixD* EV = new TMatrixD(cov_eigen.GetEigenVectors());
    m_EV.push_back(EV);
    m_MeanValues.push_back(MeanValues);
    m_SigmaValues.push_back(SigmaValues);
    m_Gauss_means.push_back(Gauss_means);
    m_Gauss_rms.push_back(Gauss_rms);

    std::vector<std::string> layer;
    const std::vector<int>& layerNr = m_RelevantLayers;

    for (unsigned int i = 0; i < layerNr.size(); i++) {
      layer.emplace_back(Form("layer%i", layerNr[i]));
    }
    layer.emplace_back("totalE");

    std::vector<TFCS1DFunction*> cumulative;
    cumulative.reserve(layer.size());

    for (unsigned int l = 0; l < layer.size(); l++) {
      file->cd(Form("%s/bin%i/%s", folder.c_str(), bin, layer[l].c_str()));

      TFCS1DFunction* fct;
      fct = (TFCS1DFunction*)gDirectory->Get("TFCS1DFunctionRegression");
      if (!fct)
        fct = (TFCS1DFunction*)gDirectory->Get("TFCS1DFunctionRegressionTF");
      if (!fct)
        fct = (TFCS1DFunction*)gDirectory->Get("TFCS1DFunctionHistogram");
      cumulative.push_back(fct);
    }

    m_cumulative.emplace_back(std::move(cumulative));
  }
  m_totalE_probability_ratio.resize(m_numberpcabins - 1, nullptr);

  return true;
}

void TFCSPCAEnergyParametrization::clean()
{
  for (unsigned int i = 0; i < m_EV.size(); i++)
    delete m_EV[i];
}

void TFCSPCAEnergyParametrization::Streamer(TBuffer& R__b)
{
  // Stream an object of class TFCSPCAEnergyParametrization

  if (R__b.IsReading()) {
    UInt_t R__s, R__c;
    Version_t R__v = R__b.ReadVersion(&R__s, &R__c);
    R__b.SetBufferOffset(R__s);

    R__b.ReadClassBuffer(TFCSPCAEnergyParametrization::Class(), this);

    if (R__v == 1) {
      set_total_energy_normalization(Ekin_nominal());
      // Check if min and max range is a factor 2 within 1 MeV tolerance
      // If yes, the nominal to log-avergage
      if (TMath::Abs(Ekin_max() - 2 * Ekin_min()) < 1) {
        set_Ekin_nominal(Ekin_max() / TMath::Sqrt(2));
      }
    }
    if (R__v <= 2) {
      m_totalE_probability_ratio.resize(m_numberpcabins - 1, nullptr);
    }
  } else {
    R__b.WriteClassBuffer(TFCSPCAEnergyParametrization::Class(), this);
  }
}
