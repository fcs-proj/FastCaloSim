/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#ifndef TFCSHitCellMappingWiggle_h
#define TFCSHitCellMappingWiggle_h

#include "FastCaloSim/TFCSHitCellMapping.h"

class TFCS1DFunction;
class TH1;

class TFCSHitCellMappingWiggle : public TFCSHitCellMapping
{
public:
  TFCSHitCellMappingWiggle(const char* name = nullptr,
                           const char* title = nullptr,
                           ICaloGeometry* geo = nullptr);
  ~TFCSHitCellMappingWiggle();

  void initialize(TFCS1DFunction* func);
  void initialize(const std::vector<const TFCS1DFunction*>& functions,
                  const std::vector<float>& bin_low_edges);

  void initialize(TH1* histogram, float xscale = 1);
  void initialize(const std::vector<const TH1*>& histograms,
                  const std::vector<float>& bin_low_edges,
                  float xscale = 1);

  inline unsigned int get_number_of_bins() const { return m_functions.size(); };

  inline double get_bin_low_edge(int bin) const { return m_bin_low_edge[bin]; };
  inline double get_bin_up_edge(int bin) const
  {
    return m_bin_low_edge[bin + 1];
  };

  inline const TFCS1DFunction* get_function(int bin) const
  {
    return m_functions[bin];
  };
  const std::vector<const TFCS1DFunction*> get_functions()
  {
    return m_functions;
  };
  const std::vector<float> get_bin_low_edges() { return m_bin_low_edge; };

  /// modify one hit position to emulate the LAr accordion shape
  /// and then fills all hits into calorimeter cells
  virtual FCSReturnCode simulate_hit(
      Hit& hit,
      TFCSSimulationState& simulstate,
      const TFCSTruthState* truth,
      const TFCSExtrapolationState* extrapol) override;

  virtual bool operator==(const TFCSParametrizationBase& ref) const override;

  void Print(Option_t* option = "") const override;

  static void unit_test(TFCSSimulationState* simulstate = nullptr,
                        TFCSTruthState* truth = nullptr,
                        TFCSExtrapolationState* extrapol = nullptr);

protected:
  bool compare(const TFCSParametrizationBase& ref) const;

private:
  //** Function for the hit-to-cell assignment accordion structure fix (wiggle)
  //**//
  //** To be moved to the conditions database at some point **//
  std::vector<const TFCS1DFunction*> m_functions = {nullptr};
  std::vector<float> m_bin_low_edge = {0, static_cast<float>(init_eta_max)};

  ClassDefOverride(TFCSHitCellMappingWiggle, 1)  // TFCSHitCellMappingWiggle
};

#endif
