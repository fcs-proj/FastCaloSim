// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project
#pragma once

#include <cstdlib>
#include <fstream>
#include <mutex>
#include <unordered_set>

// Debug-only instrumentation helper (FCS_DEBUG_EXTRAPOL_EVENT investigation).
// inline + function-local static gives a single mutex instance shared across
// every translation unit that includes this header, so prints from different
// source files (CaloGeo.cxx, FastCaloSimCaloExtrapolation.cxx) serialize
// against each other too, not just against themselves. Without this, MT's
// concurrent worker threads interleave std::cout's operator<< chains
// mid-line and corrupt the printed values.
inline std::mutex& fcsDebugPrintMutex()
{
  static std::mutex m;
  return m;
}

// FCS_DEBUG_PHI_FILTER_FILE=<path>: optional file of IDCaloBoundary_phi
// values (one per line, e.g. from an earlier FCS_DEBUG_EXTRAPOL_EVENT pass),
// used to restrict debug output to just the particles of one
// already-identified event instead of every event in the job -- cuts output
// volume by roughly the number of events without losing any per-particle
// detail for the event actually being investigated. Shared via this header
// (rather than duplicated per-file) so CaloGeo.cxx's RTree probe can filter
// on the same key as FastCaloSimCaloExtrapolation.cxx's prints: the
// extrapolateToLayers() reference-surface lookup queries CaloGeo with
// pos.phi() set to exactly IDCaloBoundary_phi(), so the same filter value
// applies to both. Doubles round-trip exactly through a
// 17-significant-digit text representation (IEEE 754), so exact-match
// lookup is safe here as long as both sides print with
// std::setprecision(17). Absent env var / empty file => filter is empty =>
// print everything (backwards compatible with the unfiltered
// instrumentation).
inline const std::unordered_set<double>& fcsDebugPhiFilter()
{
  static const std::unordered_set<double> filter = []
  {
    std::unordered_set<double> s;
    const char* path = std::getenv("FCS_DEBUG_PHI_FILTER_FILE");
    if (path) {
      std::ifstream in(path);
      double v;
      while (in >> v)
        s.insert(v);
    }
    return s;
  }();
  return filter;
}

inline bool fcsDebugShouldPrint(double phi)
{
  const auto& filter = fcsDebugPhiFilter();
  return filter.empty() || filter.count(phi) > 0;
}
