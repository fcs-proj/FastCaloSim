// Copyright (c) 2024 CERN for the benefit of the FastCaloSim project

#include "FastCaloSim/Core/TFCSEnergyBinParametrization.h"

#include "CLHEP/Random/RandFlat.h"
#include "FastCaloSim/Core/TFCSSimulationState.h"
#include "FastCaloSim/Core/TFCSTruthState.h"
#include "TFile.h"
#include "TMath.h"
#include "TVectorF.h"

//=============================================
//======= TFCSSelectEnergyBin =========
//=============================================
using namespace FastCaloSim::Core;

TFCSEnergyBinParametrization::TFCSEnergyBinParametrization(const char* name,
                                                           const char* title)
    : TFCSEnergyParametrization(name, title)
{
  m_number_of_Ekin_bins = 1;
}

bool TFCSEnergyBinParametrization::is_match_Ekin_bin(int Ekin_bin) const
{
  if (Ekin_bin >= 1 && Ekin_bin <= n_bins())
    return true;
  return false;
}

void TFCSEnergyBinParametrization::resize()
{
  for (int id : pdgid()) {
    MSG_VERBOSE("PDGid=" << id << " resize to " << n_bins() + 1);
    m_pdgid_Ebin_probability[id].resize(n_bins() + 1);
  }
  for (auto it = m_pdgid_Ebin_probability.begin();
       it != m_pdgid_Ebin_probability.end();
       ++it)
  {
    if (!is_match_pdgid(it->first)) {
      m_pdgid_Ebin_probability.erase(it);
      it = m_pdgid_Ebin_probability.begin();
    }
  }
}

void TFCSEnergyBinParametrization::set_pdgid(int id)
{
  TFCSEnergyParametrization::set_pdgid(id);
  resize();
}

void TFCSEnergyBinParametrization::set_pdgid(const std::set<int>& ids)
{
  TFCSEnergyParametrization::set_pdgid(ids);
  resize();
}

void TFCSEnergyBinParametrization::add_pdgid(int id)
{
  TFCSEnergyParametrization::add_pdgid(id);
  resize();
}

void TFCSEnergyBinParametrization::clear_pdgid()
{
  TFCSEnergyParametrization::clear_pdgid();
  resize();
}

void TFCSEnergyBinParametrization::set_pdgid_Ekin_bin_probability(
    int id, std::vector<float> prob)
{
  add_pdgid(id);
  if (prob.size() != m_pdgid_Ebin_probability[id].size()) {
    MSG_ERROR(
        "TFCSEnergyBinParametrization::set_pdgid_Ekin_bin_"
        "probability(): size of vectors does not match! in.size()="
        << prob.size() << " instance=" << m_pdgid_Ebin_probability[id].size());
    return;
  }
  float ptot = 0;
  for (int iEbin = 0; iEbin <= n_bins(); ++iEbin)
    ptot += prob[iEbin];
  float p = 0;
  for (int iEbin = 0; iEbin <= n_bins(); ++iEbin) {
    p += prob[iEbin] / ptot;
    m_pdgid_Ebin_probability[id][iEbin] = p;
  }
}

bool TFCSEnergyBinParametrization::load_pdgid_Ekin_bin_probability_from_file(
    int id, TFile* file, std::string prob_object_name)
{
  add_pdgid(id);

  file->cd();
  auto* probFromFile =
      dynamic_cast<TVectorF*>(gDirectory->Get(prob_object_name.c_str()));

  std::vector<float> prob;
  prob.reserve(m_pdgid_Ebin_probability[id].size());
  if (probFromFile == nullptr) {
    MSG_INFO(
        "TFCSEnergyBinParametrization::load_pdgid_Ekin_bin_"
        "probability_from_file(): "
        << prob_object_name << " is null. Using equal PCA probabilities.");
    prob.push_back(0.0);
    for (size_t i = 1; i < m_pdgid_Ebin_probability[id].size(); i++) {
      prob.push_back(1.0);
    }
  } else {
    auto size = static_cast<size_t>(probFromFile->GetNoElements());
    if (size != m_pdgid_Ebin_probability[id].size()) {
      MSG_ERROR(
          "TFCSEnergyBinParametrization::load_pdgid_Ekin_bin_probability_from_"
          "file(): size of prob array does not match! in.size()="
          << size << " instance=" << m_pdgid_Ebin_probability[id].size());
      return false;
    }

    const float* probArray = probFromFile->GetMatrixArray();
    for (size_t i = 0; i < size; i++) {
      prob.push_back(probArray[i]);
    }
  }

  float ptot {};
  for (int iEbin = 0; iEbin <= n_bins(); ++iEbin)
    ptot += prob[iEbin];
  float p {};
  for (int iEbin = 0; iEbin <= n_bins(); ++iEbin) {
    p += prob[iEbin] / ptot;
    m_pdgid_Ebin_probability[id][iEbin] = p;
  }

  return true;
}

void TFCSEnergyBinParametrization::Print(Option_t* option) const
{
  TString opt(option);
  bool shortprint = opt.Index("short") >= 0;
  bool longprint = msgLvl(MSG::DEBUG) || (msgLvl(MSG::INFO) && !shortprint);
  TString optprint = opt;
  optprint.ReplaceAll("short", "");
  TFCSEnergyParametrization::Print(option);
  if (longprint) {
    for (std::set<int>::iterator it = pdgid().begin(); it != pdgid().end();
         ++it)
    {
      MSG(INFO) << optprint << "  PDGID=" << *it << " : ";
      float p = 0;
      for (int iEbin = 0; iEbin <= n_bins(); ++iEbin) {
        if (iEbin > 0)
          MSG(INFO) << ", ";
        auto mapit = m_pdgid_Ebin_probability.find(*it);
        MSG(INFO) << "b" << iEbin << "="
                  << (mapit->second[iEbin] - p) / mapit->second.back() * 100
                  << "%";
        p = mapit->second[iEbin];
      }
      MSG(INFO) << END_MSG(INFO);
    }
  }
}

FCSReturnCode TFCSEnergyBinParametrization::simulate(
    TFCSSimulationState& simulstate,
    const TFCSTruthState* truth,
    const TFCSExtrapolationState* /*extrapol*/) const
{
  if (!simulstate.randomEngine()) {
    return FCSFatal;
  }

  const int truth_pdgid = truth->pdgid();
  int pdgid = -99;

  if (is_match_pdgid(truth_pdgid))
    pdgid = truth_pdgid;
  else if (is_match_pdgid(0))
    pdgid = 0;
  else {
    MSG_ERROR("TFCSEnergyBinParametrization::simulate(): cannot simulate pdgid="
              << truth_pdgid);
    return FCSFatal;
  }

  const float searchRand = CLHEP::RandFlat::shoot(simulstate.randomEngine());
  const auto& Ebin_probability = m_pdgid_Ebin_probability.at(pdgid);
  int chosenBin =
      TMath::BinarySearch(n_bins() + 1, Ebin_probability.data(), searchRand)
      + 1;
  if (chosenBin < 0) {
    MSG_WARNING(
        "TFCSEnergyBinParametrization::simulate(): chosenBin<0 "
        "(will use chosenBin=0)");
    std::string array = "";
    for (const auto& prob : Ebin_probability) {
      array += prob;
      array += " ";
    }
    MSG_WARNING(" E=" << simulstate.E() << " Ebin=" << chosenBin
                      << " rnd=" << searchRand << " array=" << array);
    chosenBin = 0;
  } else if (chosenBin > n_bins()) {
    MSG_WARNING(
        "TFCSEnergyBinParametrization::simulate(): "
        "chosenBin>n_bins() (will use chosenBin=n_bins())");
    std::string array = "";
    for (const auto& prob : Ebin_probability) {
      array += prob;
      array += " ";
    }
    MSG_WARNING(" E=" << simulstate.E() << " Ebin=" << chosenBin
                      << " rnd=" << searchRand << " array=" << array);
    chosenBin = n_bins();
  }
  simulstate.set_Ebin(chosenBin);
  MSG_DEBUG("Ebin=" << chosenBin);

  return FCSSuccess;
}
