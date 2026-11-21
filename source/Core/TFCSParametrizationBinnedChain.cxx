// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#include <format>
#include <iostream>

#include "FastCaloSim/Core/TFCSParametrizationBinnedChain.h"

#include "FastCaloSim/Core/TFCSExtrapolationState.h"
#include "FastCaloSim/Core/TFCSSimulationState.h"
#include "FastCaloSim/Core/TFCSTruthState.h"

//=============================================
//======= TFCSParametrizationBinnedChain =========
//=============================================

void TFCSParametrizationBinnedChain::push_before_first_bin(
    TFCSParametrizationBase* param)
{
  if (m_bin_start[0] == size()) {
    chain().push_back(param);
  } else {
    Chain_t::iterator const it(&chain()[m_bin_start[0]]);
    chain().insert(it, param);
  }
  for (auto& ele : m_bin_start) {
    ele++;
  }
  recalc();
}

void TFCSParametrizationBinnedChain::push_back_in_bin(
    TFCSParametrizationBase* param, std::size_t bin)
{
  if (m_bin_start.size() <= bin + 1) {
    m_bin_start.resize(bin + 2, m_bin_start.back());
    m_bin_start.shrink_to_fit();
  }
  if (m_bin_start[bin + 1] == size()) {
    chain().push_back(param);
  } else {
    Chain_t::iterator const it(&chain()[m_bin_start[bin + 1]]);
    chain().insert(it, param);
  }
  for (std::size_t ibin = bin + 1; ibin < m_bin_start.size(); ++ibin) {
    m_bin_start[ibin]++;
  }
  recalc();
}

auto TFCSParametrizationBinnedChain::get_bin(
    TFCSSimulationState& /*unused*/,
    const TFCSTruthState* /*unused*/,
    const TFCSExtrapolationState* /*unused*/) const
    -> std::optional<std::size_t>
{
  return 0;
}

auto TFCSParametrizationBinnedChain::find_bin(
    TFCSSimulationState& simulstate,
    const TFCSTruthState* truth,
    const TFCSExtrapolationState* extrapol) const -> std::optional<std::size_t>
{
  const auto bin = get_bin(simulstate, truth, extrapol);
  if (bin && *bin < get_number_of_bins()) {
    return bin;
  }
  return std::nullopt;
}

auto TFCSParametrizationBinnedChain::bin_label(std::optional<std::size_t> bin)
    -> std::string
{
  return bin ? std::to_string(*bin) : std::string("none");
}

// NOLINTNEXTLINE(readability-const-return-type): virtual, overridden
auto TFCSParametrizationBinnedChain::get_variable_text(
    TFCSSimulationState& /*unused*/,
    const TFCSTruthState* /*unused*/,
    const TFCSExtrapolationState* /*unused*/) const -> const std::string
{
  return std::string("NO VARIABLE DEFINED");
}

// NOLINTNEXTLINE(readability-const-return-type): virtual, overridden
auto TFCSParametrizationBinnedChain::get_bin_text(
    std::optional<std::size_t> bin) const -> const std::string
{
  return std::format("bin {}", bin_label(bin));
}

auto TFCSParametrizationBinnedChain::simulate(
    TFCSSimulationState& simulstate,
    const TFCSTruthState* truth,
    const TFCSExtrapolationState* extrapol) const -> FCSReturnCode
{
  return simulate_with_restarts(
      "TFCSParametrizationBinnedChain::simulate()",
      [&] -> FCSReturnCode
      {
        // Parametrizations before the bins, run for all bins
        for (std::size_t ichain = 0; ichain < m_bin_start[0]; ++ichain) {
          FCS_MSG_DEBUG("now run for all bins: " << chain()[ichain]->GetName());
          const FCSReturnCode status =
              simulate_and_retry(chain()[ichain], simulstate, truth, extrapol);
          if (status != FCSSuccess) {
            return status;
          }
        }

        // Parametrizations of the bin of this particle
        if (get_number_of_bins() > 0) {
          const auto bin = find_bin(simulstate, truth, extrapol);
          if (bin) {
            for (unsigned int ichain = m_bin_start[*bin];
                 ichain < m_bin_start[*bin + 1];
                 ++ichain)
            {
              FCS_MSG_DEBUG("for "
                            << get_variable_text(simulstate, truth, extrapol)
                            << " run " << get_bin_text(bin) << ": "
                            << chain()[ichain]->GetName());
              const FCSReturnCode status = simulate_and_retry(
                  chain()[ichain], simulstate, truth, extrapol);
              if (status != FCSSuccess) {
                return status;
              }
            }
          } else {
            FCS_MSG_WARNING(
                "for " << get_variable_text(simulstate, truth, extrapol) << ": "
                       << get_bin_text(get_bin(simulstate, truth, extrapol)));
          }
        } else {
          FCS_MSG_WARNING("no bins defined, is this intended?");
        }

        // Parametrizations after the bins, run for all bins
        for (std::size_t ichain = m_bin_start.back(); ichain < size(); ++ichain)
        {
          FCS_MSG_DEBUG("now run for all bins: " << chain()[ichain]->GetName());
          const FCSReturnCode status =
              simulate_and_retry(chain()[ichain], simulstate, truth, extrapol);
          if (status != FCSSuccess) {
            return status;
          }
        }
        return FCSSuccess;
      });
}

void TFCSParametrizationBinnedChain::Print(Option_t* option) const
{
  // The chain is printed per bin below, not as in TFCSParametrizationChain
  // NOLINTNEXTLINE(bugprone-parent-virtual-call): skip the chain printout
  TFCSParametrization::Print(option);
  const std::string opt = option != nullptr ? option : "";
  const auto [shortprint, longprint, optprint] = print_options(option);

  std::string prefix = "- ";
  for (std::size_t ichain = 0; ichain < size(); ++ichain) {
    if (ichain == 0 && ichain != m_bin_start.front()) {
      prefix = "> ";
      if (longprint) {
        FCS_MSG_INFO(optprint << prefix << "Run for all bins");
      }
    }
    for (std::size_t ibin = 0; ibin < get_number_of_bins(); ++ibin) {
      if (ichain == m_bin_start[ibin]) {
        if (ibin < get_number_of_bins() - 1) {
          if (ichain == m_bin_start[ibin + 1]) {
            continue;
          }
        }
        prefix = std::format("{:<2}", ibin);
        if (longprint) {
          FCS_MSG_INFO(optprint << prefix << "Run for " << get_bin_text(ibin));
        }
      }
    }
    if (ichain == m_bin_start.back()) {
      prefix = "< ";
      if (longprint) {
        FCS_MSG_INFO(optprint << prefix << "Run for all bins");
      }
    }
    chain()[ichain]->Print((opt + prefix).c_str());
  }
}
