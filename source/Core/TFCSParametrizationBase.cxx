// Copyright (c) 2024 CERN for the benefit of the FastCaloSim project

#include "FastCaloSim/Core/TFCSParametrizationBase.h"

#include "FastCaloSim/Geometry/CaloGeo.h"
#include "TClass.h"

//=============================================
//======= TFCSParametrizationBase =========
//=============================================

TFCSParametrizationBase::TFCSParametrizationBase(const char* name,
                                                 const char* title)
    : TNamed(name, title)
{
}

void TFCSParametrizationBase::set_geometry(CaloGeo* geo)
{
  for (unsigned int i = 0; i < size(); ++i)
    (*this)[i]->set_geometry(geo);
}

/// Result should be returned in simulstate.
/// Simulate all energies in calo layers for energy parametrizations.
/// Simulate cells for shape simulation.
FCSReturnCode TFCSParametrizationBase::simulate(
    TFCSSimulationState& /*simulstate*/,
    const TFCSTruthState* /*truth*/,
    const TFCSExtrapolationState* /*extrapol*/) const
{
  FCS_MSG_ERROR(
      "now in TFCSParametrizationBase::simulate(). This should "
      "normally not happen");
  // Force one retry to issue a printout from the chain causing the call to this
  // method
  return (FCSReturnCode)(FCSRetry + 1);
}

bool TFCSParametrizationBase::compare(const TFCSParametrizationBase& ref) const
{
  if (this == &ref) {
    FCS_MSG_DEBUG("compare(): identical instances " << this << " == " << &ref);
    return true;
  }
  return false;
}

/// If called with argument "short", only a one line summary will be printed
void TFCSParametrizationBase::Print(Option_t* option) const
{
  TString opt(option);
  bool shortprint = opt.Index("short") >= 0;
  bool longprint =
      msgLvl(FCS_MSG::DEBUG) || (msgLvl(FCS_MSG::INFO) && !shortprint);
  TString optprint = opt;
  optprint.ReplaceAll("short", "");

  if (longprint) {
    FCS_MSG_INFO(optprint << GetTitle() << " " << IsA()->GetName());
    FCS_MSG(INFO) << optprint << "  PDGID: ";
    if (is_match_all_pdgid()) {
      FCS_MSG(INFO) << "all";
    } else {
      for (std::set<int>::iterator it = pdgid().begin(); it != pdgid().end();
           ++it)
      {
        if (it != pdgid().begin())
          FCS_MSG(INFO) << ", ";
        FCS_MSG(INFO) << *it;
      }
    }
    if (is_match_all_Ekin()) {
      FCS_MSG(INFO) << " ; Ekin=all";
    } else {
      FCS_MSG(INFO) << " ; Ekin=" << Ekin_nominal() << " [" << Ekin_min()
                    << " , " << Ekin_max() << ") MeV";
    }
    if (is_match_all_eta()) {
      FCS_MSG(INFO) << " ; eta=all";
    } else {
      FCS_MSG(INFO) << " ; eta=" << eta_nominal() << " [" << eta_min() << " , "
                    << eta_max() << ")";
    }
    FCS_MSG(INFO) << END_FCS_MSG(INFO);
  } else {
    FCS_MSG_INFO(optprint << GetTitle());
  }
}

void TFCSParametrizationBase::FindDuplicates(FindDuplicateClasses_t& dupclasses)
{
  for (unsigned int i = 0; i < size(); ++i)
    if ((*this)[i]) {
      TFCSParametrizationBase* param = (*this)[i];
      FindDuplicates_t& dup = dupclasses[param->GetName()];
      // If param is already in the duplication list, skip over
      auto checkexist = dup.find(param);
      if (checkexist != dup.end()) {
        FCS_MSG_DEBUG("Found duplicate pointer for: " << param << "="
                                                      << param->GetName());
        if (checkexist->second.replace) {
          TFCSParametrizationBase* refparam = checkexist->second.replace;
          FCS_MSG_DEBUG("Found duplicate pointer: "
                        << refparam << "=" << refparam->GetName()
                        << ", duplicate is " << param << "=" << param->GetName()
                        << " index " << i << " of " << this);
          dup[refparam].mother.push_back(this);
          dup[refparam].index.push_back(i);
        }
        continue;
      }
      // Add param to duplication list
      dup[param] = Duplicate_t();
      for (auto& ref : dup) {
        TFCSParametrizationBase* refparam = ref.first;
        // skip itself, as this just added above
        if (param == refparam)
          continue;
        // skip nullptr reference
        if (refparam == nullptr)
          continue;
        // skip reference that is itself going to get replaced
        if (ref.second.replace)
          continue;
        // Check for objects with identical content
        if (*param == *refparam) {
          FCS_MSG_DEBUG("Found duplicate: "
                        << refparam << "=" << refparam->GetName()
                        << ", duplicate is " << param << "=" << param->GetName()
                        << " index " << i << " of " << this);
          dup[param].replace = refparam;
          dup[refparam].mother.push_back(this);
          dup[refparam].index.push_back(i);
          break;
        }
      }
      // Continue for child objects in param
      param->FindDuplicates(dupclasses);
    }
}

void TFCSParametrizationBase::RemoveDuplicates()
{
  FindDuplicateClasses_t dupclasses;
  FindDuplicates(dupclasses);

  std::set<TFCSParametrizationBase*> dellist;
  for (auto& dupiter : dupclasses) {
    FindDuplicates_t& dup = dupiter.second;
    for (auto onedup : dup) {
      if (onedup.second.mother.empty())
        continue;
      TFCSParametrizationBase* ref = onedup.first;
      FCS_MSG_DEBUG("Main object " << ref << "=" << ref->GetName());
      for (unsigned int i = 0; i < onedup.second.mother.size(); ++i) {
        int index = onedup.second.index[i];
        TFCSParametrizationBase* mother = onedup.second.mother[i];
        TFCSParametrizationBase* delparam = mother->operator[](index);
        unsigned int delcount = dup[delparam].mother.size();
        if (delcount == 0) {
          FCS_MSG_DEBUG("  - Delete object "
                        << delparam << "=" << delparam->GetName() << " index "
                        << index << " of " << mother << ", has " << delcount
                        << " other replacements attached. Deleting");
          mother->set_daughter(index, ref);
          dellist.insert(delparam);
        } else {
          FCS_MSG_WARNING("  - Delete object "
                          << delparam << "=" << delparam->GetName() << " index "
                          << index << " of " << mother << ", has " << delcount
                          << " other replacements attached. Skipping");
        }
      }
    }
  }

  FCS_MSG_INFO("RERUNNING DUPLICATE FINDING");
  FindDuplicateClasses_t dupclasses2;
  FindDuplicates(dupclasses2);

  std::map<std::string, int> ndel;
  for (auto* delparam : dellist) {
    FindDuplicates_t& dup2 = dupclasses2[delparam->GetName()];
    bool present = dup2.find(delparam) != dup2.end();
    if (present) {
      FCS_MSG_WARNING("- Delete object " << delparam << "="
                                         << delparam->GetName()
                                         << " still referenced somewhere!");
    } else {
      FCS_MSG_DEBUG("- Delete object " << delparam << "="
                                       << delparam->GetName());
      ++ndel[delparam->ClassName()];
      delete delparam;
    }
  }
  for (auto& del : ndel)
    FCS_MSG_INFO("Deleted " << del.second << " duplicate objects of class "
                            << del.first);
}

void TFCSParametrizationBase::RemoveNameTitle()
{
  for (unsigned int i = 0; i < size(); ++i)
    if ((*this)[i]) {
      TFCSParametrizationBase* param = (*this)[i];
      param->SetName("");
      param->SetTitle("");

      // Continue for child objects in param
      param->RemoveNameTitle();
    }
}
