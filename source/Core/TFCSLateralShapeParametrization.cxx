// Copyright (c) 2024 CERN for the benefit of the FastCaloSim project

#include "FastCaloSim/Core/TFCSLateralShapeParametrization.h"

#include <TClass.h>

//=============================================
//======= TFCSLateralShapeParametrization =========
//=============================================

TFCSLateralShapeParametrization::TFCSLateralShapeParametrization(
    const char* name, const char* title)
    : TFCSParametrization(name, title)
    , m_Ekin_bin(-1)
    , m_calosample(-1)
{
}

void TFCSLateralShapeParametrization::set_Ekin_bin(int bin)
{
  m_Ekin_bin = bin;
}

void TFCSLateralShapeParametrization::set_calosample(int cs)
{
  m_calosample = cs;
}

void TFCSLateralShapeParametrization::set_pdgid_Ekin_eta_Ekin_bin_calosample(
    const TFCSLateralShapeParametrization& ref)
{
  set_calosample(ref.calosample());
  set_Ekin_bin(ref.Ekin_bin());
  set_pdgid_Ekin_eta(ref);
}

bool TFCSLateralShapeParametrization::compare(
    const TFCSParametrizationBase& ref) const
{
  if (IsA() != ref.IsA()) {
    FCS_MSG_DEBUG("compare(): different class types "
                  << IsA()->GetName() << " != " << ref.IsA()->GetName());
    return false;
  }
  const TFCSLateralShapeParametrization& ref_typed =
      static_cast<const TFCSLateralShapeParametrization&>(ref);
  if (Ekin_bin() != ref_typed.Ekin_bin()) {
    FCS_MSG_DEBUG("compare(): different Ekin bin");
    return false;
  }
  if (calosample() != ref_typed.calosample()) {
    FCS_MSG_DEBUG("compare(): different calosample");
    return false;
  }

  return true;
}

void TFCSLateralShapeParametrization::Print(Option_t* option) const
{
  TString opt(option);
  bool shortprint = opt.Index("short") >= 0;
  bool longprint =
      msgLvl(FCS_MSG::DEBUG) || (msgLvl(FCS_MSG::INFO) && !shortprint);
  TString optprint = opt;
  optprint.ReplaceAll("short", "");
  TFCSParametrization::Print(option);
  if (longprint) {
    if (Ekin_bin() == -1)
      FCS_MSG_INFO(optprint << "  Ekin_bin=all ; calosample=" << calosample());
    else
      FCS_MSG_INFO(optprint << "  Ekin_bin=" << Ekin_bin()
                            << " ; calosample=" << calosample());
  }
}
