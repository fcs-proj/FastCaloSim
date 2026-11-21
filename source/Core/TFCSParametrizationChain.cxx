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
  return {param.Ekin_nominal(), param.Ekin_min(), param.Ekin_max()};
}

auto eta_range(const TFCSParametrizationBase& param) -> Range
{
  return {param.eta_nominal(), param.eta_min(), param.eta_max()};
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
  const auto nominal_outside = [&range]
  { return range.nominal < range.min || range.nominal > range.max; };
  for (const auto* param : chain) {
    const Range other = range_of(*param);
    if (intersect ? other.min > range.min : other.min < range.min)
      range.min = other.min;
    if (intersect ? other.max < range.max : other.max > range.max)
      range.max = other.max;
    if (nominal_outside())
      range.nominal = other.nominal;
  }
  if (nominal_outside())
    range.nominal = 0.5 * (range.min + range.max);
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
    if (intersect)
      std::ranges::set_intersection(
          pdgids, param->pdgid(), std::inserter(combined, combined.begin()));
    else
      std::ranges::set_union(
          pdgids, param->pdgid(), std::inserter(combined, combined.begin()));
    pdgids = std::move(combined);
  }
  return pdgids;
}
}  // namespace

void TFCSParametrizationChain::recalc_pdgid_intersect()
{
  set_pdgid(combined_pdgids(m_chain, true));
}

void TFCSParametrizationChain::recalc_pdgid_union()
{
  set_pdgid(combined_pdgids(m_chain, false));
}

void TFCSParametrizationChain::recalc_Ekin_intersect()
{
  set_Ekin_range(*this, combined_range(m_chain, Ekin_range, true));
}

void TFCSParametrizationChain::recalc_eta_intersect()
{
  set_eta_range(*this, combined_range(m_chain, eta_range, true));
}

void TFCSParametrizationChain::recalc_Ekin_eta_intersect()
{
  recalc_Ekin_intersect();
  recalc_eta_intersect();
}

void TFCSParametrizationChain::recalc_Ekin_union()
{
  set_Ekin_range(*this, combined_range(m_chain, Ekin_range, false));
}

void TFCSParametrizationChain::recalc_eta_union()
{
  set_eta_range(*this, combined_range(m_chain, eta_range, false));
}

void TFCSParametrizationChain::recalc_Ekin_eta_union()
{
  recalc_Ekin_union();
  recalc_eta_union();
}

void TFCSParametrizationChain::recalc()
{
  clear();
  if (m_chain.empty())
    return;

  recalc_pdgid_intersect();
  recalc_Ekin_eta_intersect();

  m_chain.shrink_to_fit();
}

auto TFCSParametrizationChain::is_match_Ekin_bin(int Ekin_bin) const -> bool
{
  for (const auto& param : m_chain)
    if (!param->is_match_Ekin_bin(Ekin_bin))
      return false;
  return true;
}

auto TFCSParametrizationChain::is_match_calosample(int calosample) const -> bool
{
  for (const auto& param : m_chain)
    if (!param->is_match_calosample(calosample))
      return false;
  return true;
}

auto TFCSParametrizationChain::simulate(
    TFCSSimulationState& simulstate,
    const TFCSTruthState* truth,
    const TFCSExtrapolationState* extrapol) const -> FCSReturnCode
{
  Int_t retry = 0;
  Int_t retry_warning = 1;

  FCSReturnCode status = FCSSuccess;
  for (int i = 0; i <= retry; i++) {
    if (i >= retry_warning)
      FCS_MSG_WARNING(
          "TFCSParametrizationChain::simulate(): Retry simulate call "
          << i << "/" << retry);
    for (const auto& param : m_chain) {
      status = simulate_and_retry(param, simulstate, truth, extrapol);

      if (status >= FCSRetry) {
        retry = status - FCSRetry;
        retry_warning = retry >> 1;
        if (retry_warning < 1)
          retry_warning = 1;
        break;
      }
      if (status == FCSFatal)
        return FCSFatal;
    }

    if (status == FCSSuccess)
      break;
  }

  if (status != FCSSuccess) {
    FCS_MSG_FATAL(
        "TFCSParametrizationChain::simulate(): Simulate call failed after "
        << retry << " retries");
    return FCSFatal;
  }

  return FCSSuccess;
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

void TFCSParametrizationChain::Streamer(TBuffer& R__b)
{
  // Stream an object of class TFCSParametrizationChain.

  UInt_t R__s, R__c;
  TDirectory* dir = nullptr;

  if (R__b.IsReading()) {
    Version_t R__v = R__b.ReadVersion(&R__s, &R__c);
    if (R__v == 1) {
      R__b.SetBufferOffset(R__s);
      R__b.ReadClassBuffer(TFCSParametrizationChain::Class(), this);
    } else {
      TFCSParametrization::Streamer(R__b);

      TObject* parent = R__b.GetParent();
      if (R__b.GetParent()) {
        if (parent->InheritsFrom(TDirectory::Class())) {
          dir = static_cast<TDirectory*>(parent);
        }
      }

      TFCSParametrizationChain::Chain_t& R__stl = m_chain;
      R__stl.clear();
      TClass* R__tcl1 = TFCSParametrizationBase::Class();
      if (R__tcl1 == nullptr) {
        Error("m_chain streamer",
              "Missing the TClass object for class TFCSParametrizationBase *!");
        return;
      }
      int R__i, R__n;
      R__b >> R__n;
      R__stl.reserve(R__n);
      for (R__i = 0; R__i < R__n; R__i++) {
        std::unique_ptr<TFCSParametrizationBase> R__t;
        R__t.reset((TFCSParametrizationBase*)R__b.ReadObjectAny(R__tcl1));
        if (R__t != nullptr) {
          if (R__t->InheritsFrom(TFCSParametrizationPlaceholder::Class())) {
            std::unique_ptr<TFCSParametrizationBase> new_R__t = nullptr;

            if (dir) {
              new_R__t.reset(
                  (TFCSParametrizationBase*)dir->Get(R__t->GetName()));
            }

            if (new_R__t) {
              R__t = std::move(new_R__t);
            } else {
              Error("TFCSParametrizationChain::Streamer",
                    "Found placeholder object in the parametrization chain, "
                    "but could not read the real object from the file!");
            }
          }
        }
        if (R__t) {
          R__t->fixHists();
        } else {
          Error("TFCSParametrizationChain::Streamer",
                "Failed to read element %d of %d in the parametrization chain",
                R__i,
                R__n);
        }
        R__stl.emplace_back(R__t.release());
      }

      R__b.CheckByteCount(R__s, R__c, TFCSParametrizationChain::IsA());
    }
  } else {
    R__c = R__b.WriteVersion(TFCSParametrizationChain::IsA(), kTRUE);
    TFCSParametrization::Streamer(R__b);

    if (SplitChainObjects()) {
      TObject* parent = R__b.GetParent();
      if (R__b.GetParent()) {
        if (parent->InheritsFrom(TDirectory::Class())) {
          dir = static_cast<TDirectory*>(parent);
        }
      }
    }

    TFCSParametrizationChain::Chain_t& R__stl = m_chain;
    int R__n = int(R__stl.size());
    R__b << R__n;
    if (R__n) {
      TFCSParametrizationChain::Chain_t::iterator R__k;
      for (R__k = R__stl.begin(); R__k != R__stl.end(); ++R__k) {
        TFCSParametrizationBase* R__t = *R__k;  // Ownership stays with m_chain
        if (dir && R__t != nullptr) {
          dir->WriteTObject(R__t);
          TFCSParametrizationPlaceholder tmp(
              R__t->GetName(), TString("Placeholder for: ") + R__t->GetTitle());
          R__b.WriteObject(&tmp, false);  // tell R__b object with same address
                                          // are actually different
        } else {
          R__b.WriteObject(R__t);
        }
      }
    }
    R__b.SetByteCount(R__c, kTRUE);
  }
}
