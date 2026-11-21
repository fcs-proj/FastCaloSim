// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#include <algorithm>
#include <iterator>
#include <set>
#include <utility>

#include "FastCaloSim/Core/TFCSParametrizationChain.h"

#include "FastCaloSim/Core/TFCSExtrapolationState.h"
#include "FastCaloSim/Core/TFCSParametrizationPlaceholder.h"
#include "FastCaloSim/Core/TFCSSimulationState.h"
#include "FastCaloSim/Core/TFCSTruthState.h"
#include "TBuffer.h"
#include "TDirectory.h"

namespace
{
// Kinetic energy or eta range that a parametrization covers
struct Range
{
  double nominal;
  double min;
  double max;
};

auto Ekin_range(const TFCSParametrizationBase& param) -> Range
{
  return {.nominal = param.Ekin_nominal(),
          .min = param.Ekin_min(),
          .max = param.Ekin_max()};
}

auto eta_range(const TFCSParametrizationBase& param) -> Range
{
  return {.nominal = param.eta_nominal(),
          .min = param.eta_min(),
          .max = param.eta_max()};
}

// Intersection or union of the ranges of all parametrizations in a chain. The
// nominal value of the first parametrization is kept as long as it lies inside
// the combined range, otherwise the nominal value of the parametrization just
// added is taken. If it is still outside at the end, the center is taken.
auto combined_range(const TFCSParametrizationChain::Chain_t& chain,
                    Range (*range_of)(const TFCSParametrizationBase&),
                    bool intersect) -> Range
{
  Range range = range_of(*chain[0]);
  const auto nominal_outside = [&range] -> bool
  { return range.nominal < range.min || range.nominal > range.max; };
  for (const auto* param : chain) {
    const Range other = range_of(*param);
    if (intersect ? other.min > range.min : other.min < range.min) {
      range.min = other.min;
    }
    if (intersect ? other.max < range.max : other.max > range.max) {
      range.max = other.max;
    }
    if (nominal_outside()) {
      range.nominal = other.nominal;
    }
  }
  if (nominal_outside()) {
    range.nominal = 0.5 * (range.min + range.max);
  }
  return range;
}

void set_Ekin_range(TFCSParametrization& param, const Range& range)
{
  param.set_Ekin_nominal(range.nominal);
  param.set_Ekin_min(range.min);
  param.set_Ekin_max(range.max);
}

void set_eta_range(TFCSParametrization& param, const Range& range)
{
  param.set_eta_nominal(range.nominal);
  param.set_eta_min(range.min);
  param.set_eta_max(range.max);
}

// Intersection or union of the particle IDs of all parametrizations in a chain
auto combined_pdgids(const TFCSParametrizationChain::Chain_t& chain,
                     bool intersect) -> std::set<int>
{
  std::set<int> pdgids = chain[0]->pdgid();
  for (const auto* param : chain) {
    std::set<int> combined;
    if (intersect) {
      std::ranges::set_intersection(
          pdgids, param->pdgid(), std::inserter(combined, combined.begin()));
    } else {
      std::ranges::set_union(
          pdgids, param->pdgid(), std::inserter(combined, combined.begin()));
    }
    pdgids = std::move(combined);
  }
  return pdgids;
}
}  // namespace

void TFCSParametrizationChain::recalc_pdgid_intersect()
{
  set_pdgid(combined_pdgids(m_chain, /*intersect=*/true));
}

void TFCSParametrizationChain::recalc_pdgid_union()
{
  set_pdgid(combined_pdgids(m_chain, /*intersect=*/false));
}

void TFCSParametrizationChain::recalc_Ekin_intersect()
{
  set_Ekin_range(*this,
                 combined_range(m_chain, Ekin_range, /*intersect=*/true));
}

void TFCSParametrizationChain::recalc_eta_intersect()
{
  set_eta_range(*this, combined_range(m_chain, eta_range, /*intersect=*/true));
}

void TFCSParametrizationChain::recalc_Ekin_eta_intersect()
{
  recalc_Ekin_intersect();
  recalc_eta_intersect();
}

void TFCSParametrizationChain::recalc_Ekin_union()
{
  set_Ekin_range(*this,
                 combined_range(m_chain, Ekin_range, /*intersect=*/false));
}

void TFCSParametrizationChain::recalc_eta_union()
{
  set_eta_range(*this, combined_range(m_chain, eta_range, /*intersect=*/false));
}

void TFCSParametrizationChain::recalc_Ekin_eta_union()
{
  recalc_Ekin_union();
  recalc_eta_union();
}

void TFCSParametrizationChain::recalc()
{
  clear();
  if (m_chain.empty()) {
    return;
  }

  recalc_pdgid_intersect();
  recalc_Ekin_eta_intersect();

  m_chain.shrink_to_fit();
}

auto TFCSParametrizationChain::is_match_Ekin_bin(std::size_t Ekin_bin) const
    -> bool
{
  return std::ranges::all_of(
      m_chain,
      [Ekin_bin](const TFCSParametrizationBase* param) -> bool
      { return param->is_match_Ekin_bin(Ekin_bin); });
}

auto TFCSParametrizationChain::is_match_calosample(std::size_t calosample) const
    -> bool
{
  return std::ranges::all_of(
      m_chain,
      [calosample](const TFCSParametrizationBase* param) -> bool
      { return param->is_match_calosample(calosample); });
}

auto TFCSParametrizationChain::simulate(
    TFCSSimulationState& simulstate,
    const TFCSTruthState* truth,
    const TFCSExtrapolationState* extrapol) const -> FCSReturnCode
{
  return simulate_with_restarts("TFCSParametrizationChain::simulate()",
                                [&] -> FCSReturnCode
                                {
                                  for (const auto& param : m_chain) {
                                    const FCSReturnCode status =
                                        simulate_and_retry(
                                            param, simulstate, truth, extrapol);
                                    if (status != FCSSuccess) {
                                      return status;
                                    }
                                  }
                                  return FCSSuccess;
                                });
}

void TFCSParametrizationChain::Print(Option_t* option) const
{
  TFCSParametrization::Print(option);
  const std::string opt = option != nullptr ? option : "";

  char count = 'A';
  for (const auto& param : m_chain) {
    param->Print((opt + count + ' ').c_str());
    count++;
  }
}

void TFCSParametrizationChain::Streamer(TBuffer& buffer)
{
  // Stream an object of class TFCSParametrizationChain.

  UInt_t start = 0;
  UInt_t count = 0;
  TDirectory* dir = nullptr;

  if (buffer.IsReading()) {
    Version_t const version = buffer.ReadVersion(&start, &count);
    if (version == 1) {
      buffer.SetBufferOffset(static_cast<Int_t>(start));
      buffer.ReadClassBuffer(TFCSParametrizationChain::Class(), this);
    } else {
      TFCSParametrization::Streamer(buffer);

      TObject* parent = buffer.GetParent();
      if (buffer.GetParent() != nullptr) {
        if (parent->InheritsFrom(TDirectory::Class())) {
          dir = static_cast<TDirectory*>(parent);
        }
      }

      TFCSParametrizationChain::Chain_t& params = m_chain;
      params.clear();
      TClass const* base_class = TFCSParametrizationBase::Class();
      if (base_class == nullptr) {
        Error("m_chain streamer",
              "Missing the TClass object for class TFCSParametrizationBase *!");
        return;
      }
      int i = 0;
      int n = 0;
      buffer >> n;
      params.reserve(static_cast<std::size_t>(n));
      for (i = 0; i < n; i++) {
        std::unique_ptr<TFCSParametrizationBase> param;
        param.reset(static_cast<TFCSParametrizationBase*>(
            buffer.ReadObjectAny(base_class)));
        if (param != nullptr) {
          if (param->InheritsFrom(TFCSParametrizationPlaceholder::Class())) {
            std::unique_ptr<TFCSParametrizationBase> new_param = nullptr;

            if (dir != nullptr) {
              new_param.reset(static_cast<TFCSParametrizationBase*>(
                  dir->Get(param->GetName())));
            }

            if (new_param) {
              param = std::move(new_param);
            } else {
              Error("TFCSParametrizationChain::Streamer",
                    "Found placeholder object in the parametrization chain, "
                    "but could not read the real object from the file!");
            }
          }
        }
        if (param) {
          param->fixHists();
        } else {
          Error("TFCSParametrizationChain::Streamer",
                "Failed to read element %d of %d in the parametrization chain",
                i,
                n);
        }
        params.emplace_back(param.release());
      }

      buffer.CheckByteCount(start, count, TFCSParametrizationChain::IsA());
    }
  } else {
    count = buffer.WriteVersion(TFCSParametrizationChain::IsA(), kTRUE);
    TFCSParametrization::Streamer(buffer);

    if (SplitChainObjects()) {
      TObject* parent = buffer.GetParent();
      if (buffer.GetParent() != nullptr) {
        if (parent->InheritsFrom(TDirectory::Class())) {
          dir = static_cast<TDirectory*>(parent);
        }
      }
    }

    TFCSParametrizationChain::Chain_t& params = m_chain;
    int const n = static_cast<int>(params.size());
    buffer << n;
    if (n != 0) {
      TFCSParametrizationChain::Chain_t::iterator it;
      for (it = params.begin(); it != params.end(); ++it) {
        TFCSParametrizationBase const* param =
            *it;  // Ownership stays with m_chain
        if (dir != nullptr && param != nullptr) {
          dir->WriteTObject(param);
          TFCSParametrizationPlaceholder const tmp(
              param->GetName(),
              TString("Placeholder for: ") + param->GetTitle());
          buffer.WriteObject(&tmp, false);  // tell buffer object with same
                                            // address are actually different
        } else {
          buffer.WriteObject(param);
        }
      }
    }
    buffer.SetByteCount(count, kTRUE);
  }
}
