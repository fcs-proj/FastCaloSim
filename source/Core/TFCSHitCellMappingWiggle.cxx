#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-equal"

// Copyright (c) 2024 CERN for the benefit of the FastCaloSim project

#include "FastCaloSim/Core/TFCSHitCellMappingWiggle.h"

#include <TClass.h>

#include "CLHEP/Random/RandFlat.h"
#include "FastCaloSim/Core/TFCS1DFunctionInt32Histogram.h"
#include "FastCaloSim/Core/TFCSExtrapolationState.h"
#include "FastCaloSim/Core/TFCSSimulationState.h"
#include "FastCaloSim/Core/TFCSTruthState.h"
#include "TH1.h"
#include "TVector2.h"

//=============================================
//======= TFCSHitCellMappingWiggle =========
//=============================================

TFCSHitCellMappingWiggle::TFCSHitCellMappingWiggle(const char* name,
                                                   const char* title,
                                                   CaloGeo* geo)
    : TFCSHitCellMapping(name, title, geo)
{
}

TFCSHitCellMappingWiggle::~TFCSHitCellMappingWiggle()
{
  for (const auto* function : m_functions)
    delete function;
}

void TFCSHitCellMappingWiggle::initialize(TFCS1DFunction* func)
{
  if (!func)
    return;
  for (const auto* function : m_functions)
    if (function)
      delete function;

  m_functions.resize(1);
  m_functions[0] = func;

  m_bin_low_edge.resize(2);
  m_bin_low_edge[0] = 0;
  m_bin_low_edge[1] = init_eta_max;
}

void TFCSHitCellMappingWiggle::initialize(
    const std::vector<const TFCS1DFunction*>& functions,
    const std::vector<float>& bin_low_edges)
{
  if (functions.size() + 1 != bin_low_edges.size()) {
    ATH_MSG_ERROR("Using " << functions.size() << " functions needs "
                           << functions.size() + 1 << " bin low edges, but got "
                           << bin_low_edges.size() << "bins");
    return;
  }
  for (const auto* function : m_functions)
    if (function)
      delete function;
  m_functions = functions;
  m_bin_low_edge = bin_low_edges;
}

void TFCSHitCellMappingWiggle::initialize(TH1* histogram, float xscale)
{
  if (!histogram)
    return;
  TFCS1DFunctionInt32Histogram* func =
      new TFCS1DFunctionInt32Histogram(histogram);
  if (xscale != 1) {
    for (auto& ele : func->get_HistoBordersx())
      ele *= xscale;
  }
  initialize(func);
}

void TFCSHitCellMappingWiggle::initialize(
    const std::vector<const TH1*>& histograms,
    const std::vector<float>& bin_low_edges,
    float xscale)
{
  if (histograms.size() + 1 != bin_low_edges.size()) {
    ATH_MSG_ERROR("Using " << histograms.size() << " histograms needs "
                           << histograms.size() + 1 << " bins, but got "
                           << bin_low_edges.size() << "bins");
    return;
  }
  std::vector<const TFCS1DFunction*> functions(histograms.size());
  for (unsigned int i = 0; i < histograms.size(); ++i) {
    if (histograms[i]) {
      TFCS1DFunctionInt32Histogram* func =
          new TFCS1DFunctionInt32Histogram(histograms[i]);
      if (xscale != 1) {
        for (auto& ele : func->get_HistoBordersx())
          ele *= xscale;
      }
      functions[i] = func;
    } else {
      functions[i] = nullptr;
    }
  }

  initialize(functions, bin_low_edges);
}

FCSReturnCode TFCSHitCellMappingWiggle::simulate_hit(
    Hit& hit,
    TFCSSimulationState& simulstate,
    const TFCSTruthState* truth,
    const TFCSExtrapolationState* extrapol)
{
  if (!simulstate.randomEngine()) {
    return FCSFatal;
  }

  float eta = fabs(hit.eta());
  if (eta < m_bin_low_edge[0] || eta >= m_bin_low_edge[get_number_of_bins()]) {
    return TFCSHitCellMapping::simulate_hit(hit, simulstate, truth, extrapol);
  }

  auto it = std::upper_bound(m_bin_low_edge.begin(), m_bin_low_edge.end(), eta);
  int bin = std::distance(m_bin_low_edge.begin(), it) - 1;

  const TFCS1DFunction* func = get_function(bin);
  if (func) {
    double rnd = CLHEP::RandFlat::shoot(simulstate.randomEngine());

    double wiggle = func->rnd_to_fct(rnd);

    ATH_MSG_DEBUG("HIT: E=" << hit.E() << " cs=" << calosample()
                            << " eta=" << hit.eta() << " phi=" << hit.phi()
                            << " wiggle=" << wiggle << " bin=" << bin << " ["
                            << get_bin_low_edge(bin) << ","
                            << get_bin_up_edge(bin) << "] func=" << func);

    double hit_phi_shifted = hit.phi() + wiggle;
    hit.set_phi_y(TVector2::Phi_mpi_pi(hit_phi_shifted));
  }

  return TFCSHitCellMapping::simulate_hit(hit, simulstate, truth, extrapol);
}

bool TFCSHitCellMappingWiggle::operator==(
    const TFCSParametrizationBase& ref) const
{
  if (TFCSParametrizationBase::compare(ref))
    return true;
  if (!TFCSParametrization::compare(ref))
    return false;
  if (!TFCSLateralShapeParametrization::compare(ref))
    return false;
  if (!TFCSHitCellMappingWiggle::compare(ref))
    return false;

  return true;
}

void TFCSHitCellMappingWiggle::Print(Option_t* option) const
{
  TFCSHitCellMapping::Print(option);
  TString opt(option);
  bool shortprint = opt.Index("short") >= 0;
  bool longprint = msgLvl(MSG::DEBUG) || (msgLvl(MSG::INFO) && !shortprint);
  TString optprint = opt;
  optprint.ReplaceAll("short", "");

  if (longprint) {
    ATH_MSG(INFO) << optprint << "  " << get_number_of_bins()
                  << " functions : ";
    for (unsigned int i = 0; i < get_number_of_bins(); ++i)
      ATH_MSG(INFO) << get_bin_low_edge(i) << " < (" << get_function(i)
                    << ") < ";
    ATH_MSG(INFO) << get_bin_up_edge(get_number_of_bins() - 1) << END_MSG(INFO);
  }
}

bool TFCSHitCellMappingWiggle::compare(const TFCSParametrizationBase& ref) const
{
  if (IsA() != ref.IsA()) {
    ATH_MSG_DEBUG("compare(): different class types "
                  << IsA()->GetName() << " != " << ref.IsA()->GetName());
    return false;
  }
  const TFCSHitCellMappingWiggle& ref_typed =
      static_cast<const TFCSHitCellMappingWiggle&>(ref);

  if (m_bin_low_edge != ref_typed.m_bin_low_edge) {
    ATH_MSG_DEBUG("operator==(): different bin edges");
    return false;
  }

  for (unsigned int i = 0; i < get_number_of_bins(); ++i) {
    const TFCS1DFunction* f1 = get_function(i);
    const TFCS1DFunction* f2 = ref_typed.get_function(i);
    if (!f1 && !f2)
      continue;
    if ((f1 && !f2) || (!f1 && f2)) {
      ATH_MSG_DEBUG(
          "compare(): different only one function pointer is nullptr");
      return false;
    }
    if (f1->IsA() != f2->IsA()) {
      ATH_MSG_DEBUG("compare(): different class types for function "
                    << i << ": " << f1->IsA()->GetName()
                    << " != " << f2->IsA()->GetName());
      return false;
    }
    if (!(*f1 == *f2)) {
      ATH_MSG_DEBUG("compare(): difference in functions " << i);
      return false;
    }
  }

  return true;
}

#pragma GCC diagnostic pop
