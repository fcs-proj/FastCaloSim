// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#ifndef TFCSHitCellMappingWiggle_h
#define TFCSHitCellMappingWiggle_h

#include "FastCaloSim/Core/TFCSHitCellMapping.h"

class TFCS1DFunction;
class TH1;

class TFCSHitCellMappingWiggle : public TFCSHitCellMapping
{
public:
  TFCSHitCellMappingWiggle(const char* name = nullptr,
                           const char* title = nullptr,
                           CaloGeo* geo = nullptr);
  ~TFCSHitCellMappingWiggle() override;

  void initialize(TFCS1DFunction* func);
  void initialize(const std::vector<const TFCS1DFunction*>& functions,
                  const std::vector<float>& bin_low_edges);

  void initialize(TH1* histogram, float xscale = 1);
  void initialize(const std::vector<const TH1*>& histograms,
                  const std::vector<float>& bin_low_edges,
                  float xscale = 1);

  void clear();

  inline auto get_number_of_bins() const -> unsigned int
  {
    return m_functions.size();
  };

  inline auto get_bin_low_edge(int bin) const -> double
  {
    return m_bin_low_edge[bin];
  };
  inline auto get_bin_up_edge(int bin) const -> double
  {
    return m_bin_low_edge[bin + 1];
  };

  inline auto get_function(int bin) const -> const TFCS1DFunction*
  {
    return m_functions[bin];
  };
  auto get_functions() const -> const std::vector<const TFCS1DFunction*>&
  {
    return m_functions;
  };
  auto get_bin_low_edges() const -> const std::vector<float>&
  {
    return m_bin_low_edge;
  };

  /// modify one hit position to emulate the LAr accordion shape
  /// and then fills all hits into calorimeter cells
  auto simulate_hit(Hit& hit,
                    TFCSSimulationState& simulstate,
                    const TFCSTruthState* truth,
                    const TFCSExtrapolationState* extrapol)
      -> FCSReturnCode override;

  auto operator==(const TFCSParametrizationBase& ref) const -> bool override;

  void Print(Option_t* option = "") const override;

protected:
  auto compare(const TFCSParametrizationBase& ref) const -> bool;

private:
  //** Function for the hit-to-cell assignment accordion structure fix (wiggle)
  //**//
  //** To be moved to the conditions database at some point **//
  std::vector<const TFCS1DFunction*> m_functions = {nullptr};
  std::vector<float> m_bin_low_edge = {0, static_cast<float>(init_eta_max)};

  ClassDefOverride(TFCSHitCellMappingWiggle, 1)  // TFCSHitCellMappingWiggle
};

#endif
