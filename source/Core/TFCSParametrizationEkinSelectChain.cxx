// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#include <cmath>
#include <cstddef>
#include <format>
#include <iostream>

#include "FastCaloSim/Core/TFCSParametrizationEkinSelectChain.h"

#include "CLHEP/Random/RandFlat.h"
#include "FastCaloSim/Core/TFCSExtrapolationState.h"
#include "FastCaloSim/Core/TFCSSimulationState.h"
#include "FastCaloSim/Core/TFCSTruthState.h"

//=============================================
//======= TFCSParametrizationEkinSelectChain =========
//=============================================

void TFCSParametrizationEkinSelectChain::recalc()
{
  clear();
  if (size() == 0) {
    return;
  }

  recalc_pdgid_intersect();
  recalc_Ekin_union();
  recalc_eta_intersect();

  chain().shrink_to_fit();
}

void TFCSParametrizationEkinSelectChain::push_back_in_bin(
    TFCSParametrizationBase* param)
{
  // The bin edges are stored in single precision
  push_back_in_bin(param,
                   static_cast<float>(param->Ekin_min()),
                   static_cast<float>(param->Ekin_max()));
}

auto TFCSParametrizationEkinSelectChain::get_bin(
    TFCSSimulationState& simulstate,
    const TFCSTruthState* truth,
    const TFCSExtrapolationState*) const -> std::optional<std::size_t>
{
  if (!simulstate.randomEngine()) {
    return std::nullopt;
  }

  const double Ekin = truth->Ekin();
  // The bin edges are single precision, so is the binning variable
  const auto found = val_to_bin(static_cast<float>(Ekin));

  if (!DoRandomInterpolation()) {
    return found;
  }
  // Drawn also when no bin is found, to keep the random number sequence
  const double rnd = CLHEP::RandFlat::shoot(simulstate.randomEngine());
  if (!found) {
    return found;
  }
  std::size_t bin = *found;

  // First parametrization in a bin, nullptr if the bin is empty
  const auto first_param_in =
      [this](std::size_t ibin) -> TFCSParametrizationBase*
  {
    if (m_bin_start[ibin + 1] == m_bin_start[ibin]) {
      return nullptr;
    }
    return chain()[m_bin_start[ibin]];
  };

  // if no parametrizations for this bin, return
  TFCSParametrizationBase* first_in_bin = first_param_in(bin);
  if (!first_in_bin) {
    return bin;
  }

  if (Ekin < first_in_bin->Ekin_nominal()) {
    if (bin == 0) {
      return bin;
    }
    const std::size_t prevbin = bin - 1;
    // if no parametrizations for previous bin, return
    TFCSParametrizationBase* first_in_prevbin = first_param_in(prevbin);
    if (!first_in_prevbin) {
      return bin;
    }

    const double logEkin = std::log(Ekin);
    const double logEkin_nominal = std::log(first_in_bin->Ekin_nominal());
    const double logEkin_previous = std::log(first_in_prevbin->Ekin_nominal());
    const double numerator = logEkin - logEkin_previous;
    const double denominator = logEkin_nominal - logEkin_previous;
    if (denominator <= 0) {
      return bin;
    }

    if (numerator / denominator < rnd) {
      bin = prevbin;
    }
    FCS_MSG_DEBUG("logEkin="
                  << logEkin << " logEkin_previous=" << logEkin_previous
                  << " logEkin_nominal=" << logEkin_nominal
                  << " (rnd=" << 1 - rnd
                  << " < p(previous)=" << (1 - numerator / denominator)
                  << ")? => orgbin=" << prevbin + 1 << " selbin=" << bin);
  } else {
    if (bin + 1 == get_number_of_bins()) {
      return bin;
    }
    const std::size_t nextbin = bin + 1;
    // if no parametrizations for next bin, return
    TFCSParametrizationBase* first_in_nextbin = first_param_in(nextbin);
    if (!first_in_nextbin) {
      return bin;
    }

    const double logEkin = std::log(Ekin);
    const double logEkin_nominal = std::log(first_in_bin->Ekin_nominal());
    const double logEkin_next = std::log(first_in_nextbin->Ekin_nominal());
    const double numerator = logEkin - logEkin_nominal;
    const double denominator = logEkin_next - logEkin_nominal;
    if (denominator <= 0) {
      return bin;
    }

    if (rnd < numerator / denominator) {
      bin = nextbin;
    }
    FCS_MSG_DEBUG("logEkin="
                  << logEkin << " logEkin_nominal=" << logEkin_nominal
                  << " logEkin_next=" << logEkin_next << " (rnd=" << rnd
                  << " < p(next)=" << numerator / denominator
                  << ")? => orgbin=" << nextbin - 1 << " selbin=" << bin);
  }

  return bin;
}

auto TFCSParametrizationEkinSelectChain::get_variable_text(
    TFCSSimulationState&,
    const TFCSTruthState* truth,
    const TFCSExtrapolationState*) const -> const std::string
{
  return std::format("Ekin={:1.1f}", truth->Ekin());
}

auto TFCSParametrizationEkinSelectChain::get_bin_text(
    std::optional<std::size_t> bin) const -> const std::string
{
  if (!bin || *bin >= get_number_of_bins()) {
    return std::format("bin={} not in [{:1.1f}<=Ekin<{:1.1f})",
                       bin_label(bin),
                       m_bin_low_edge[0],
                       m_bin_low_edge[get_number_of_bins()]);
  }
  if (DoRandomInterpolation()) {
    return std::format("bin={}, {:1.1f}<=Ekin(+random)<{:1.1f}",
                       *bin,
                       get_bin_low_edge(*bin),
                       get_bin_up_edge(*bin));
  }
  return std::format("bin={}, {:1.1f}<=Ekin<{:1.1f}",
                     *bin,
                     get_bin_low_edge(*bin),
                     get_bin_up_edge(*bin));
}
