// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#include <algorithm>
#include <cstddef>
#include <iostream>
#include <iterator>

#include "FastCaloSim/Core/TFCSParametrizationFloatSelectChain.h"

#include "FastCaloSim/Definitions/FloatCompare.h"
#include "FastCaloSim/Definitions/Index.h"

//=============================================
//======= TFCSParametrizationFloatSelectChain =========
//=============================================

auto TFCSParametrizationFloatSelectChain::push_back_in_bin(
    TFCSParametrizationBase* param, float low, float up)
    -> std::optional<std::size_t>
{
  if (up < low) {
    // can't handle wrong order of bounds
    FCS_MSG_ERROR("Cannot add " << param->GetName() << ": range [" << low << ","
                                << up << ") not well defined");
    return std::nullopt;
  }
  if (get_number_of_bins() == 0) {
    // special case of adding the first bin
    push_back_in_bin(param, 0);
    m_bin_low_edge[0] = low;
    m_bin_low_edge[1] = up;
    return 0;
  }
  const auto ilow = val_to_bin(low);
  const auto iup = val_to_bin(up);
  if (!ilow && !iup
      && FastCaloSim::exactly_equal(m_bin_low_edge[get_number_of_bins()], low))
  {
    // special case of adding bin at the end of existing bins
    const std::size_t endbin = get_number_of_bins();
    push_back_in_bin(param, endbin);
    m_bin_low_edge[endbin] = low;
    m_bin_low_edge[endbin + 1] = up;
    return endbin;
  }

  if (!ilow && !iup) {
    // can't handle disjunct ranges
    FCS_MSG_ERROR("Cannot add " << param->GetName() << ": range [" << low << ","
                                << up << ") which is outside existing range ["
                                << m_bin_low_edge[0] << ","
                                << m_bin_low_edge[get_number_of_bins()] << ")");
    return std::nullopt;
  }

  if (ilow && iup && *iup == *ilow + 1
      && FastCaloSim::exactly_equal(m_bin_low_edge[*ilow], low)
      && FastCaloSim::exactly_equal(m_bin_low_edge[*iup], up))
  {
    // Case of adding to an existing bin
    push_back_in_bin(param, *ilow);
    return ilow;
  }

  if (!ilow && iup == 0 && FastCaloSim::exactly_equal(m_bin_low_edge[0], up)) {
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
  return std::nullopt;
}

void TFCSParametrizationFloatSelectChain::push_back_in_bin(
    TFCSParametrizationBase* param, std::size_t bin)
{
  TFCSParametrizationBinnedChain::push_back_in_bin(param, bin);
  if (m_bin_low_edge.size() < m_bin_start.size()) {
    m_bin_low_edge.resize(m_bin_start.size(), m_bin_low_edge.back());
    m_bin_low_edge.shrink_to_fit();
  }
}

auto TFCSParametrizationFloatSelectChain::val_to_bin(float val) const
    -> std::optional<std::size_t>
{
  if (val < m_bin_low_edge[0]) {
    FCS_MSG_VERBOSE("val_to_bin(" << val << ")=none: " << val << " < "
                                  << m_bin_low_edge[0]);
    return std::nullopt;
  }
  if (val >= m_bin_low_edge[get_number_of_bins()]) {
    FCS_MSG_VERBOSE("val_to_bin(" << val << ")=none: " << val << " >= "
                                  << m_bin_low_edge[get_number_of_bins()]);
    return std::nullopt;
  }

  // val is within the edges, so the edge above it is not the first one
  const auto it = std::ranges::upper_bound(m_bin_low_edge, val);
  const auto bin = FastCaloSim::index_of(m_bin_low_edge, it) - 1;
  FCS_MSG_VERBOSE("val_to_bin(" << val << ")=" << bin);
  return bin;
}
