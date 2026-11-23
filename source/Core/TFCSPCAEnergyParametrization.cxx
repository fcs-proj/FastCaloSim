// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <format>
#include <memory>
#include <numbers>
#include <utility>
#include <vector>

#include "FastCaloSim/Core/TFCSPCAEnergyParametrization.h"

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGaussZiggurat.h"
#include "FastCaloSim/Core/GaussZiggurat.h"
#include "FastCaloSim/Core/TFCSExtrapolationState.h"
#include "FastCaloSim/Core/TFCSSimulationState.h"
#include "FastCaloSim/Geometry/CaloGeo.h"
#include "TClass.h"
#include "TFile.h"
#include "TH1.h"
#include "TKey.h"
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
}

auto TFCSPCAEnergyParametrization::is_match_Ekin_bin(std::size_t Ekin_bin) const
    -> bool
{
  return Ekin_bin >= 1 && Ekin_bin <= n_bins();
}

auto TFCSPCAEnergyParametrization::is_match_calosample(
    std::size_t calosample) const -> bool
{
  return std::ranges::any_of(m_RelevantLayers,
                             [calosample](int layer) -> bool
                             { return std::cmp_equal(layer, calosample); });
}

void TFCSPCAEnergyParametrization::Print(Option_t* option) const
{
  const auto [shortprint, longprint, optprint] = print_options(option);
  TFCSEnergyParametrization::Print(option);

  if (longprint) {
    FCS_MSG(INFO) << optprint << "  #bins=" << m_numberpcabins
                  << ", Enorm=" << m_total_energy_normalization << ", layers=";
    for (std::size_t i = 0; i < m_RelevantLayers.size(); ++i) {
      if (i > 0) {
        FCS_MSG(INFO) << ", ";
      }
      FCS_MSG(INFO) << m_RelevantLayers[i];
    }
    FCS_MSG(INFO) << END_FCS_MSG(INFO);
  }
}

namespace
{
auto interpolate_get_y(TH1* hist, double x) -> double
{
  double m = 0;
  double n = 0;
  double x1 = 1;
  double x2 = 0;
  double y1 = 0;
  double y2 = 0;
  if (x <= hist->GetBinCenter(1)) {
    return hist->GetBinContent(1);
  }
  if (x >= hist->GetBinCenter(hist->GetNbinsX())) {
    return hist->GetBinContent(hist->GetNbinsX());
  }

  int const bin = hist->FindBin(x);
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
}  // namespace

void TFCSPCAEnergyParametrization::set_totalE_probability_ratio(
    std::size_t Ekin_bin, TH1* hist)
{
  // The Ekin bins start at 1
  if (Ekin_bin < 1 || Ekin_bin > m_totalE_probability_ratio.size()) {
    return;
  }
  m_totalE_probability_ratio[Ekin_bin - 1] = hist;
}

auto TFCSPCAEnergyParametrization::get_totalE_probability_ratio(
    std::size_t Ekin_bin) const -> TH1*
{
  // The Ekin bins start at 1
  if (Ekin_bin < 1 || Ekin_bin > m_totalE_probability_ratio.size()) {
    return nullptr;
  }
  return m_totalE_probability_ratio[Ekin_bin - 1];
}

auto TFCSPCAEnergyParametrization::simulate(
    TFCSSimulationState& simulstate,
    const TFCSTruthState* /*truth*/,
    const TFCSExtrapolationState* /*extrapol*/) const -> FCSReturnCode
{
  if (simulstate.randomEngine() == nullptr) {
    return FCSFatal;
  }

  const auto pcabin = simulstate.Ebin();
  if (!pcabin) {
    FCS_MSG_ERROR("No energy bin chosen before the PCA energy simulation");
    return FCSFatal;
  }

  if (*pcabin == 0) {
    simulstate.set_E(0);
    for (std::size_t s = 0; s < m_geo->n_layers(); ++s) {
      simulstate.set_E(s, 0.0);
      simulstate.set_Efrac(s, 0.0);
    }
  } else {
    const std::size_t ibin = *pcabin - 1;
    TMatrixD* EV = m_EV[ibin];
    TVectorD* MeanValues = m_MeanValues[ibin];
    TVectorD* SigmaValues = m_SigmaValues[ibin];
    TVectorD* Gauss_means = m_Gauss_means[ibin];
    TVectorD* Gauss_rms = m_Gauss_rms[ibin];
    const std::vector<TFCS1DFunction*>& cumulative = m_cumulative[ibin];

    const std::vector<int>& layerNr = m_RelevantLayers;

    auto const* vals_gauss_means =
        static_cast<double*>(Gauss_means->GetMatrixArray());
    double const* vals_gauss_rms = Gauss_rms->GetMatrixArray();

    std::vector<double> output_data(layerNr.size() + 1);
    std::vector<double> input_data(layerNr.size() + 1);

    for (std::size_t l = 0; l <= layerNr.size(); l++) {
      double const mean = vals_gauss_means[l];
      double const rms = vals_gauss_rms[l];
      double const gauszz =
          FastCaloSim::gauss_ziggurat(simulstate.randomEngine(), mean, rms);
      input_data[l] = gauszz;
    }

    // The relevant layers plus the total energy
    const int n_variables = static_cast<int>(layerNr.size()) + 1;
    P2X(SigmaValues,
        MeanValues,
        EV,
        n_variables,
        input_data.data(),
        output_data.data(),
        n_variables);

    std::vector<double> simdata(layerNr.size() + 1);
    double sum_fraction = 0.0;
    for (std::size_t l = 0; l <= layerNr.size(); l++) {
      double const simdata_uniform =
          (std::erf(output_data[l] / std::numbers::sqrt2) + 1) / 2.0;

      simdata[l] = cumulative[l]->rnd_to_fct(simdata_uniform);

      if (l != layerNr.size()) {  // sum up the fractions, but not the totalE
        sum_fraction += simdata[l];
      }
    }

    double scalefactor = 1.0 / sum_fraction;
    if (do_rescale == 0) {
      scalefactor = 1.0;
    }

    for (std::size_t l = 0; l < layerNr.size(); l++) {
      simdata[l] *= scalefactor;
    }

    // Apply hit-and-miss reweighting of total energy response
    // This needs to run before simulstate is modified, otherwise a clean retry
    // is not possible
    TH1* h_totalE_ratio = get_totalE_probability_ratio(*pcabin);
    if (h_totalE_ratio != nullptr) {
      double pass_probability = 0;
      const double Etot = simdata[layerNr.size()];
      if (Etot > h_totalE_ratio->GetXaxis()->GetXmin()
          && Etot < h_totalE_ratio->GetXaxis()->GetXmax())
      {
        pass_probability = interpolate_get_y(h_totalE_ratio, Etot);
      }
      const double random = CLHEP::RandFlat::shoot(simulstate.randomEngine());
      if (random > pass_probability) {
        // Assuming an extreme h_totalE_ratio histogram that causes a retry
        // in 50% of the cases, 20 retries cause an accidental WARNING for
        // every 2^10=1024 simulated particles and an accidental FATAL for
        // every 2^20, which is safe even for the largest productions
        return FCSReturnCode::retry(20);
      }
    }

    double const total_energy = simdata[layerNr.size()] * simulstate.E()
        / static_cast<double>(m_total_energy_normalization);
    simulstate.set_E(total_energy);
    FCS_MSG_DEBUG("set E to total_energy=" << total_energy);

    for (std::size_t s = 0; s < m_geo->n_layers(); ++s) {
      double energyfrac = 0.0;
      for (std::size_t l = 0; l < layerNr.size(); ++l) {
        if (std::cmp_equal(layerNr[l], s)) {
          energyfrac = simdata[l];
        }
      }
      simulstate.set_Efrac(s, energyfrac);
      simulstate.set_E(s, energyfrac * total_energy);
      simulstate.set_SF(scalefactor);
    }
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
      x[i] += p[j] * gSigmaValues[i] * (gEigenVectors[i * gNVariables + j]);
    }
  }
}

auto TFCSPCAEnergyParametrization::loadInputs(TFile* file) -> bool
{
  return loadInputs(file, "");
}

auto TFCSPCAEnergyParametrization::loadInputs(TFile* file,
                                              const std::string& folder) -> bool
{
  bool load_ok = true;

  const std::string x = folder.empty() ? "bin" : folder + "/bin";
  while (std::unique_ptr<TObject>(file->Get(
      std::format("{}{}/pca/RelevantLayers", x, m_numberpcabins).c_str())))
  {
    m_numberpcabins++;
  }
  m_numberpcabins -= 1;

  file->cd((x + "1/pca").c_str());
  IntArray const* RelevantLayers =
      static_cast<IntArray*>(gDirectory->Get("RelevantLayers"));
  if (RelevantLayers == nullptr) {
    FCS_MSG_ERROR(
        "TFCSPCAEnergyParametrization::m_RelevantLayers in first "
        "pcabin is null!");
    load_ok = false;
  }

  if (!load_ok) {
    return false;
  }

  m_RelevantLayers.reserve(static_cast<std::size_t>(RelevantLayers->GetSize()));
  for (int i = 0; i < RelevantLayers->GetSize(); i++) {
    m_RelevantLayers.push_back(static_cast<int>(RelevantLayers->GetAt(i)));
  }

  for (int bin = 1; bin <= m_numberpcabins; bin++) {
    file->cd(std::format("{}{}/pca", x, bin).c_str());

    TMatrixDSym const* symCov =
        static_cast<TMatrixDSym*>(gDirectory->Get("symCov"));
    // NOLINTNEXTLINE(misc-const-correctness): stored as non-const pointer
    TVectorD* MeanValues =
        static_cast<TVectorD*>(gDirectory->Get("MeanValues"));
    // NOLINTNEXTLINE(misc-const-correctness): stored as non-const pointer
    TVectorD* SigmaValues =
        static_cast<TVectorD*>(gDirectory->Get("SigmaValues"));
    // NOLINTNEXTLINE(misc-const-correctness): stored as non-const pointer
    TVectorD* Gauss_means =
        static_cast<TVectorD*>(gDirectory->Get("Gauss_means"));
    // NOLINTNEXTLINE(misc-const-correctness): stored as non-const pointer
    TVectorD* Gauss_rms = static_cast<TVectorD*>(gDirectory->Get("Gauss_rms"));

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

    if (!load_ok) {
      return false;
    }

    TMatrixDSymEigen const cov_eigen(*symCov);
    auto* EV = new TMatrixD(cov_eigen.GetEigenVectors());
    m_EV.push_back(EV);
    m_MeanValues.push_back(MeanValues);
    m_SigmaValues.push_back(SigmaValues);
    m_Gauss_means.push_back(Gauss_means);
    m_Gauss_rms.push_back(Gauss_rms);

    std::vector<std::string> layer;
    const std::vector<int>& layerNr = m_RelevantLayers;

    layer.reserve(layerNr.size());
    for (int i : layerNr) {
      layer.emplace_back(std::format("layer{}", i));
    }
    layer.emplace_back("totalE");

    std::vector<TFCS1DFunction*> cumulative;
    cumulative.reserve(layer.size());

    for (auto& l : layer) {
      const std::string path = std::format("{}/bin{}/{}", folder, bin, l);
      // Without the folder, gDirectory would stay in the previous one
      if (!file->cd(path.c_str())) {
        FCS_MSG_ERROR("TFCSPCAEnergyParametrization: folder " << path
                                                              << " not found");
        return false;
      }

      auto* fct = static_cast<TFCS1DFunction*>(
          gDirectory->Get("TFCS1DFunctionRegression"));
      if (fct == nullptr) {
        fct = static_cast<TFCS1DFunction*>(
            gDirectory->Get("TFCS1DFunctionRegressionTF"));
      }
      if (fct == nullptr) {
        fct = static_cast<TFCS1DFunction*>(
            gDirectory->Get("TFCS1DFunctionHistogram"));
      }
      // simulate() samples every PCA bin and layer from its function
      if (fct == nullptr) {
        FCS_MSG_ERROR("TFCSPCAEnergyParametrization: no TFCS1DFunction in "
                      << path);
        return false;
      }
      cumulative.push_back(fct);
    }

    m_cumulative.emplace_back(std::move(cumulative));
  }
  m_totalE_probability_ratio.resize(
      static_cast<std::size_t>(m_numberpcabins - 1), nullptr);

  return true;
}

void TFCSPCAEnergyParametrization::clean()
{
  for (auto& i : m_EV) {
    delete i;
  }
}

void TFCSPCAEnergyParametrization::Streamer(TBuffer& buffer)
{
  // Stream an object of class TFCSPCAEnergyParametrization

  if (buffer.IsReading()) {
    UInt_t start = 0;
    UInt_t count = 0;
    Version_t const version = buffer.ReadVersion(&start, &count);
    buffer.SetBufferOffset(static_cast<Int_t>(start));

    buffer.ReadClassBuffer(TFCSPCAEnergyParametrization::Class(), this);

    if (version == 1) {
      set_total_energy_normalization(static_cast<float>(Ekin_nominal()));
      // Check if min and max range is a factor 2 within 1 MeV tolerance
      // If yes, the nominal to log-avergage
      if (std::abs(Ekin_max() - 2 * Ekin_min()) < 1) {
        set_Ekin_nominal(Ekin_max() / std::numbers::sqrt2);
      }
    }
    if (version <= 2) {
      m_totalE_probability_ratio.resize(
          static_cast<std::size_t>(m_numberpcabins - 1), nullptr);
    }
  } else {
    buffer.WriteClassBuffer(TFCSPCAEnergyParametrization::Class(), this);
  }
}
