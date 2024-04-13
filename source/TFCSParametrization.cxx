#pragma GCC diagnostic push 
#pragma GCC diagnostic ignored "-Wfloat-equal"

/*
  Copyright (C) 2002-2024 CERN for the benefit of the ATLAS collaboration
*/

#include "FastCaloSim/TFCSParametrization.h"

#include <TClass.h>

//=============================================
//======= TFCSParametrization =========
//=============================================

TFCSParametrization::TFCSParametrization(const char* name, const char* title)
    : TFCSParametrizationBase(name, title)
    , m_Ekin_nominal(init_Ekin_nominal)
    , m_Ekin_min(init_Ekin_min)
    , m_Ekin_max(init_Ekin_max)
    , m_eta_nominal(init_eta_nominal)
    , m_eta_min(init_eta_min)
    , m_eta_max(init_eta_max)
{
}

void TFCSParametrization::clear()
{
  m_pdgid.clear();
  m_Ekin_nominal = init_Ekin_nominal;
  m_Ekin_min = init_Ekin_min;
  m_Ekin_max = init_Ekin_max;
  m_eta_nominal = init_eta_nominal;
  m_eta_min = init_eta_min;
  m_eta_max = init_eta_max;
}

void TFCSParametrization::set_pdgid(int id)
{
  m_pdgid.clear();
  m_pdgid.insert(id);
}

void TFCSParametrization::set_pdgid(const std::set<int>& ids)
{
  m_pdgid = ids;
}

void TFCSParametrization::add_pdgid(int id)
{
  m_pdgid.insert(id);
}

void TFCSParametrization::clear_pdgid()
{
  m_pdgid.clear();
}

void TFCSParametrization::set_Ekin_nominal(double nominal)
{
  m_Ekin_nominal = nominal;
}

void TFCSParametrization::set_Ekin_min(double min)
{
  m_Ekin_min = min;
}

void TFCSParametrization::set_Ekin_max(double max)
{
  m_Ekin_max = max;
}

void TFCSParametrization::set_eta_nominal(double nominal)
{
  m_eta_nominal = nominal;
}

void TFCSParametrization::set_eta_min(double min)
{
  m_eta_min = min;
}

void TFCSParametrization::set_eta_max(double max)
{
  m_eta_max = max;
}

void TFCSParametrization::set_Ekin(const TFCSParametrizationBase& ref)
{
  set_Ekin_nominal(ref.Ekin_nominal());
  set_Ekin_min(ref.Ekin_min());
  set_Ekin_max(ref.Ekin_max());
}

void TFCSParametrization::set_eta(const TFCSParametrizationBase& ref)
{
  set_eta_nominal(ref.eta_nominal());
  set_eta_min(ref.eta_min());
  set_eta_max(ref.eta_max());
}

void TFCSParametrization::set_Ekin_eta(const TFCSParametrizationBase& ref)
{
  set_Ekin(ref);
  set_eta(ref);
}

void TFCSParametrization::set_pdgid_Ekin_eta(const TFCSParametrizationBase& ref)
{
  set_Ekin_eta(ref);
  set_pdgid(ref.pdgid());
}

bool TFCSParametrization::compare(const TFCSParametrizationBase& ref) const
{
  if (IsA() != ref.IsA()) {
    ATH_MSG_DEBUG("compare(): different class types "
                  << IsA()->GetName() << " != " << ref.IsA()->GetName());
    return false;
  }
  std::string name(GetName());
  if (name != ref.GetName()) {
    ATH_MSG_DEBUG("compare(): different names " << GetName()
                                                << " != " << ref.GetName());
    return false;
  }
  std::string title(GetTitle());
  if (title != ref.GetTitle()) {
    ATH_MSG_DEBUG("compare(): different titles " << GetTitle()
                                                 << " != " << ref.GetTitle());
    return false;
  }
  if (is_match_all_pdgid() != ref.is_match_all_pdgid()
      || pdgid() != ref.pdgid())
  {
    ATH_MSG_DEBUG("compare(): different pdgids");
    return false;
  }
  if (Ekin_nominal() != ref.Ekin_nominal() || Ekin_min() != ref.Ekin_min()
      || Ekin_max() != ref.Ekin_max())
  {
    ATH_MSG_DEBUG("compare(): different Ekin range");
    return false;
  }
  if (eta_nominal() != ref.eta_nominal() || eta_min() != ref.eta_min()
      || eta_max() != ref.eta_max())
  {
    ATH_MSG_DEBUG("compare(): different eta range");
    return false;
  }

  return true;
}

#pragma GCC diagnostic pop
