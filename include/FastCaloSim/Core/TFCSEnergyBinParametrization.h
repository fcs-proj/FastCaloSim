// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#pragma once

#include <map>
#include <vector>

#include "FastCaloSim/Core/TFCSEnergyParametrization.h"
#include "TFile.h"

class TFCSEnergyBinParametrization : public TFCSEnergyParametrization
{
public:
  TFCSEnergyBinParametrization(const char* name = nullptr,
                               const char* title = nullptr);

  void set_pdgid(int id) override;
  void set_pdgid(const std::set<int>& ids) override;
  void add_pdgid(int id) override;
  void clear_pdgid() override;

  auto n_bins() const -> int override { return m_number_of_Ekin_bins; };

  /// current convention is to start Ekin_bin counting at 1, to be updated to
  /// start counting with 0
  void set_number_of_Ekin_bins(int n_Ekin_bin)
  {
    m_number_of_Ekin_bins = n_Ekin_bin;
    resize();
  };

  /// set_pdgid_Ekin_bin_probability(): format of probability should be a vector
  /// of float with probability values, e.g. for equal probability of 5 Ekin
  /// bins specify a vector with {0,1,1,1,1,1} or {0,0.2,0.2,0.2,0.2,0.2} This
  /// would give a probability of 0 for Ekin_bin=0 and a probability of 20% for
  /// 1<=Ekin_bin<=5 the function will normalize probabilities automatically, if
  /// the sum of values is not 1 current convention is to start Ekin_bin
  /// counting at 1, to be updated to start counting with 0
  virtual void set_pdgid_Ekin_bin_probability(int id, std::vector<float> prob);
  virtual auto load_pdgid_Ekin_bin_probability_from_file(
      int id, TFile* file, std::string prob_object_name) -> bool;

  auto simulate(TFCSSimulationState& simulstate,
                const TFCSTruthState* truth,
                const TFCSExtrapolationState* extrapol) const
      -> FCSReturnCode override;

  auto is_match_Ekin_bin(int Ekin_bin) const -> bool override;

  void Print(Option_t* option = "") const override;

private:
  int m_number_of_Ekin_bins;
  std::map<int, std::vector<float>> m_pdgid_Ebin_probability;

  void resize();

  ClassDefOverride(TFCSEnergyBinParametrization,
                   1)  // TFCSEnergyBinParametrization
};
