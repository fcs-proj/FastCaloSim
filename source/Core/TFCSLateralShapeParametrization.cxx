// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#include "FastCaloSim/Core/TFCSLateralShapeParametrization.h"

#include <TClass.h>

//=============================================
//======= TFCSLateralShapeParametrization =========
//=============================================

TFCSLateralShapeParametrization::TFCSLateralShapeParametrization(
    const char* name, const char* title)
    : TFCSParametrization(name, title)

{
}

void TFCSLateralShapeParametrization::set_Ekin_bin(std::size_t bin)
{
  m_Ekin_bin = static_cast<int>(bin);
}

void TFCSLateralShapeParametrization::set_calosample(std::size_t cs)
{
  m_calosample = static_cast<int>(cs);
}

void TFCSLateralShapeParametrization::set_pdgid_Ekin_eta_Ekin_bin_calosample(
    const TFCSLateralShapeParametrization& ref)
{
  m_calosample = ref.m_calosample;
  m_Ekin_bin = ref.m_Ekin_bin;
  set_pdgid_Ekin_eta(ref);
}

auto TFCSLateralShapeParametrization::compare(
    const TFCSParametrizationBase& ref) const -> bool
{
  if (IsA() != ref.IsA()) {
    FCS_MSG_DEBUG("compare(): different class types "
                  << IsA()->GetName() << " != " << ref.IsA()->GetName());
    return false;
  }
  const TFCSLateralShapeParametrization& ref_typed =
      static_cast<const TFCSLateralShapeParametrization&>(ref);
  if (m_Ekin_bin != ref_typed.m_Ekin_bin) {
    FCS_MSG_DEBUG("compare(): different Ekin bin");
    return false;
  }
  if (m_calosample != ref_typed.m_calosample) {
    FCS_MSG_DEBUG("compare(): different calosample");
    return false;
  }

  return true;
}

void TFCSLateralShapeParametrization::Print(Option_t* option) const
{
  const auto [shortprint, longprint, optprint] = print_options(option);
  TFCSParametrization::Print(option);
  if (longprint) {
    if (!Ekin_bin()) {
      FCS_MSG_INFO(optprint << "  Ekin_bin=all ; calosample=" << m_calosample);
    } else {
      FCS_MSG_INFO(optprint << "  Ekin_bin=" << m_Ekin_bin
                            << " ; calosample=" << m_calosample);
    }
  }
}
