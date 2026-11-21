#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-equal"
// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#include <algorithm>
#include <cstddef>
#include <iostream>
#include <iterator>

#include "FastCaloSim/Core/TFCSParametrizationFloatSelectChain.h"

//=============================================
//======= TFCSParametrizationFloatSelectChain =========
//=============================================

auto TFCSParametrizationFloatSelectChain::push_back_in_bin(
    TFCSParametrizationBase* param, float low, float up) -> int
{
  if (up < low) {
    // can't handle wrong order of bounds
    FCS_MSG_ERROR("Cannot add " << param->GetName() << ": range [" << low << ","
                                << up << ") not well defined");
    return -1;
  }
  if (get_number_of_bins() == 0) {
    // special case of adding the first bin
    push_back_in_bin(param, 0);
    m_bin_low_edge[0] = low;
    m_bin_low_edge[1] = up;
    return 0;
  }
  // Bin edge for a bin index that is known to be valid
  const auto edge = [this](int bin) -> float&
  { return m_bin_low_edge[static_cast<std::size_t>(bin)]; };

  int ilow = val_to_bin(low);
  int iup = val_to_bin(up);
  if (ilow < 0 && iup < 0 && m_bin_low_edge[get_number_of_bins()] == low) {
    // special case of adding bin at the end of existing bins
    const unsigned int endbin = get_number_of_bins();
    push_back_in_bin(param, endbin);
    m_bin_low_edge[endbin] = low;
    m_bin_low_edge[endbin + 1] = up;
    return static_cast<int>(endbin);
  }

  if (ilow < 0 && iup < 0) {
    // can't handle disjunct ranges
    FCS_MSG_ERROR("Cannot add " << param->GetName() << ": range [" << low << ","
                                << up << ") which is outside existing range ["
                                << m_bin_low_edge[0] << ","
                                << m_bin_low_edge[get_number_of_bins()] << ")");
    return -1;
  }

  if (iup >= 0 && ilow >= 0 && iup - ilow == 1 && edge(ilow) == low
      && edge(iup) == up)
  {
    // Case of adding to an existing bin
    push_back_in_bin(param, static_cast<unsigned int>(ilow));
    return ilow;
  }

  if (ilow < 0 && iup == 0 && edge(iup) == up) {
    // Case of adding a new first bin before existing bins
    const std::size_t oldsize = m_bin_start.size();
    m_bin_start.resize(oldsize + 1, m_bin_start.back());
    m_bin_start.shrink_to_fit();
    m_bin_low_edge.resize(oldsize + 1, m_bin_low_edge.back());
    m_bin_low_edge.shrink_to_fit();
    for (std::size_t i = oldsize; i > 0; --i) {
      m_bin_start[i] = m_bin_start[i - 1];
      m_bin_low_edge[i] = m_bin_low_edge[i - 1];
    }
    m_bin_low_edge[0] = low;
    m_bin_low_edge[1] = up;
    push_back_in_bin(param, 0);
    return 0;
  }

  FCS_MSG_ERROR("Cannot add "
                << param->GetName() << ": range [" << low << "," << up
                << ") covers more than one bin in existing range ["
                << m_bin_low_edge[0] << ","
                << m_bin_low_edge[get_number_of_bins()]
                << ") or splits an existing bin");
  return -1;
}

void TFCSParametrizationFloatSelectChain::push_back_in_bin(
    TFCSParametrizationBase* param, unsigned int bin)
{
  TFCSParametrizationBinnedChain::push_back_in_bin(param, bin);
  if (m_bin_low_edge.size() < m_bin_start.size()) {
    m_bin_low_edge.resize(m_bin_start.size(), m_bin_low_edge.back());
    m_bin_low_edge.shrink_to_fit();
  }
}

auto TFCSParametrizationFloatSelectChain::val_to_bin(float val) const -> int
{
  if (val < m_bin_low_edge[0]) {
    FCS_MSG_VERBOSE("val_to_bin(" << val << ")=-1: " << val << " < "
                                  << m_bin_low_edge[0]);
    return -1;
  }
  if (val >= m_bin_low_edge[get_number_of_bins()]) {
    FCS_MSG_VERBOSE("val_to_bin(" << val << ")=-1: " << val << " >= "
                                  << m_bin_low_edge[get_number_of_bins()]);
    return -1;
  }

  auto it = std::upper_bound(m_bin_low_edge.begin(), m_bin_low_edge.end(), val);
  const auto dist =
      static_cast<int>(std::distance(m_bin_low_edge.begin(), it) - 1);
  FCS_MSG_VERBOSE("val_to_bin(" << val << ")=" << dist);
  return dist;
}

#pragma GCC diagnostic pop
