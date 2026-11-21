// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

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
  for (std::size_t i = 0; i < size(); ++i) {
    if (TFCSParametrizationBase* param = (*this)[i]) {
      param->set_geometry(geo);
    }
  }
}

/// Result should be returned in simulstate.
/// Simulate all energies in calo layers for energy parametrizations.
/// Simulate cells for shape simulation.
auto TFCSParametrizationBase::simulate(
    TFCSSimulationState& /*simulstate*/,
    const TFCSTruthState* /*truth*/,
    const TFCSExtrapolationState* /*extrapol*/) const -> FCSReturnCode
{
  FCS_MSG_ERROR(
      "now in TFCSParametrizationBase::simulate(). This should "
      "normally not happen");
  // Force one retry to issue a printout from the chain causing the call to this
  // method
  return static_cast<FCSReturnCode>(FCSRetry + 1);
}

auto TFCSParametrizationBase::compare(const TFCSParametrizationBase& ref) const
    -> bool
{
  if (this == &ref) {
    FCS_MSG_DEBUG("compare(): identical instances " << this << " == " << &ref);
    return true;
  }
  return false;
}

auto TFCSParametrizationBase::print_options(Option_t* option) const
    -> PrintOptions
{
  std::string optprint = option != nullptr ? option : "";
  const bool shortprint = optprint.contains("short");
  for (auto pos = optprint.find("short"); pos != std::string::npos;
       pos = optprint.find("short", pos))
  {
    optprint.erase(pos, 5);
  }
  const bool longprint = msgLvl(FastCaloSim::MSG::DEBUG)
      || (msgLvl(FastCaloSim::MSG::INFO) && !shortprint);
  return {shortprint, longprint, std::move(optprint)};
}

/// If called with argument "short", only a one line summary will be printed
void TFCSParametrizationBase::Print(Option_t* option) const
{
  const auto [shortprint, longprint, optprint] = print_options(option);

  if (longprint) {
    FCS_MSG_INFO(optprint << GetTitle() << " " << IsA()->GetName());
    FCS_MSG(INFO) << optprint << "  PDGID: ";
    if (is_match_all_pdgid()) {
      FCS_MSG(INFO) << "all";
    } else {
      for (auto it = pdgid().begin(); it != pdgid().end(); ++it) {
        if (it != pdgid().begin()) {
          FCS_MSG(INFO) << ", ";
        }
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
  for (std::size_t i = 0; i < size(); ++i) {
    if (TFCSParametrizationBase* param = (*this)[i]) {
      FindDuplicates_t& dup = dupclasses[param->GetName()];
      // If param is already in the duplication list, skip over
      auto checkexist = dup.find(param);
      if (checkexist != dup.end()) {
        FCS_MSG_WARNING(
            " [TFCSParametrizationBase::FindDuplicates] "
            "DUPLICATE POINTER DETECTED");
        FCS_MSG_WARNING(" - Pointer           : " << param);
        FCS_MSG_WARNING(" - Name              : " << param->GetName());
        FCS_MSG_WARNING(" - Class             : " << param->ClassName());
        FCS_MSG_WARNING(" - Occurs in parent  : " << this << " ("
                                                  << this->ClassName() << ")"
                                                  << " index=" << i);
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
        if (param == refparam) {
          continue;
        }
        // skip nullptr reference
        if (refparam == nullptr) {
          continue;
        }
        // skip reference that is itself going to get replaced
        if (ref.second.replace) {
          continue;
        }
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
}

void TFCSParametrizationBase::RemoveDuplicates()
{
  FindDuplicateClasses_t dupclasses;
  FindDuplicates(dupclasses);

  std::set<TFCSParametrizationBase*> dellist;
  for (auto& dupiter : dupclasses) {
    FindDuplicates_t& dup = dupiter.second;
    for (auto onedup : dup) {
      if (onedup.second.mother.empty()) {
        continue;
      }
      TFCSParametrizationBase* ref = onedup.first;
      FCS_MSG_DEBUG("Main object " << ref << "=" << ref->GetName());
      for (std::size_t i = 0; i < onedup.second.mother.size(); ++i) {
        const std::size_t index = onedup.second.index[i];
        TFCSParametrizationBase* mother = onedup.second.mother[i];
        TFCSParametrizationBase* delparam = mother->operator[](index);
        const std::size_t delcount = dup[delparam].mother.size();
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
    FindDuplicates_t const& dup2 = dupclasses2[delparam->GetName()];
    bool const present = dup2.contains(delparam);
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
  for (auto& del : ndel) {
    FCS_MSG_INFO("Deleted " << del.second << " duplicate objects of class "
                            << del.first);
  }
}

void TFCSParametrizationBase::RemoveNameTitle()
{
  for (std::size_t i = 0; i < size(); ++i) {
    if (TFCSParametrizationBase* param = (*this)[i]) {
      param->SetName("");
      param->SetTitle("");

      // Continue for child objects in param
      param->RemoveNameTitle();
    }
  }
}
