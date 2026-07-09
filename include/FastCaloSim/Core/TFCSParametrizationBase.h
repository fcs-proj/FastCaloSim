// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#pragma once

#include <map>
#include <mutex>
#include <set>

#include <FastCaloSim/FastCaloSim_export.h>

#include "FastCaloSim/Core/MLogging.h"

class CaloGeo;
class TFCSSimulationState;
class TFCSTruthState;
class TFCSExtrapolationState;

/** Base class for all FastCaloSim parametrizations
Functionality in derivde classes is  provided through the simulate method. The
simulate method takes a TFCSTruthState and a TFCSExtrapolationState object as
input and provides output in a TFCSSimulationState. Parametrizations contain
information on the pdgid, range in Ekin and range in eta of particles to which
they can be applied. Several basic types of parametrization exists:
- classes derived from TFCSEnergyParametrization simulate energy information
which is written into TFCSSimulationState
- classes derived from TFCSLateralShapeParametrization simulate cell level
information for specific calorimeter layers and bins "Ebin" in the energy
parametrization
- classes derived from TFCSParametrizationChain call other parametrization.
Depending on the derived class, these other parametrization are only called
under special conditions
- a special case of TFCSLateralShapeParametrization is
TFCSLateralShapeParametrizationHitBase for hit level shape simulation through
the simulate_hit method. Hit level simulation is controlled through the special
chain TFCSLateralShapeParametrizationHitChain.
*/

/// Return codes for the simulate function
enum FCSReturnCode
{
  FCSFatal = 0,
  FCSSuccess = 1,
  FCSRetry = 2
};

#define FCS_RETRY_COUNT 3

class FASTCALOSIM_EXPORT TFCSParametrizationBase
    : public TNamed
    , public FastCaloSim::MLogging
{
public:
  TFCSParametrizationBase(const char* name = nullptr,
                          const char* title = nullptr);

  /// Status bit for FCS needs
  enum FCSStatusBits
  {
    kMatchAllPDGID = BIT(
        14)  ///< Set this bit in the TObject bit field if valid for all PDGID
  };

  virtual auto is_match_pdgid(int /*id*/) const -> bool
  {
    return TestBit(kMatchAllPDGID);
  };
  virtual auto is_match_Ekin(float /*Ekin*/) const -> bool { return false; };
  virtual auto is_match_eta(float /*eta*/) const -> bool { return false; };

  virtual auto is_match_Ekin_bin(int /*Ekin_bin*/) const -> bool
  {
    return false;
  };
  virtual auto is_match_calosample(int /*calosample*/) const -> bool
  {
    return false;
  };

  virtual auto is_match_all_pdgid() const -> bool
  {
    return TestBit(kMatchAllPDGID);
  };
  virtual auto is_match_all_Ekin() const -> bool { return false; };
  virtual auto is_match_all_eta() const -> bool { return false; };
  virtual auto is_match_all_Ekin_bin() const -> bool { return false; };
  virtual auto is_match_all_calosample() const -> bool { return false; };

  virtual auto pdgid() const -> const std::set<int>&
  {
    static const std::set<int> empty;
    return empty;
  };
  virtual auto Ekin_nominal() const -> double { return init_Ekin_nominal; };
  virtual auto Ekin_min() const -> double { return init_Ekin_min; };
  virtual auto Ekin_max() const -> double { return init_Ekin_max; };
  virtual auto eta_nominal() const -> double { return init_eta_nominal; };
  virtual auto eta_min() const -> double { return init_eta_min; };
  virtual auto eta_max() const -> double { return init_eta_max; };

  virtual void set_match_all_pdgid() { SetBit(kMatchAllPDGID); };
  virtual void reset_match_all_pdgid() { ResetBit(kMatchAllPDGID); };

  /// Method to set the geometry access pointer. Loops over daughter objects if
  /// present
  virtual void set_geometry(CaloGeo* geo);

  /// Some derived classes have daughter instances of TFCSParametrizationBase
  /// objects The size() and operator[] methods give general access to these
  /// daughters
  virtual auto size() const -> unsigned int { return 0; };

  /// Some derived classes have daughter instances of TFCSParametrizationBase
  /// objects The size() and operator[] methods give general access to these
  /// daughters
  virtual auto operator[](unsigned int /*ind*/) const
      -> const TFCSParametrizationBase*
  {
    return nullptr;
  };

  /// Some derived classes have daughter instances of TFCSParametrizationBase
  /// objects The size() and operator[] methods give general access to these
  /// daughters
  virtual auto operator[](unsigned int /*ind*/) -> TFCSParametrizationBase*
  {
    return nullptr;
  };

  /// Some derived classes have daughter instances of TFCSParametrizationBase
  /// objects The set_daughter method allows to change these daughters - expert
  /// use only! The original element at this position is not deleted
  virtual void set_daughter(unsigned int /*ind*/,
                            TFCSParametrizationBase* /*param*/) {};

  /// The == operator compares the content of instances.
  /// The implementation in the base class only returns true for a comparison
  /// with itself
  virtual auto operator==(const TFCSParametrizationBase& ref) const -> bool
  {
    return compare(ref);
  };

  /// Method in all derived classes to do some simulation
  virtual auto simulate(TFCSSimulationState& simulstate,
                        const TFCSTruthState* truth,
                        const TFCSExtrapolationState* extrapol) const
      -> FCSReturnCode;

  /// Method in all derived classes to delete objects stored in the simulstate
  /// AuxInfo
  virtual void CleanAuxInfo(TFCSSimulationState& /*simulstate*/) const {};

  /// Print object information.
  void Print(Option_t* option = "") const override;

  struct Duplicate_t
  {
    TFCSParametrizationBase* replace = nullptr;
    std::vector<TFCSParametrizationBase*> mother;
    std::vector<unsigned int> index;
  };
  using FindDuplicates_t = std::map<TFCSParametrizationBase*, Duplicate_t>;
  using FindDuplicateClasses_t = std::map<std::string, FindDuplicates_t>;
  void FindDuplicates(FindDuplicateClasses_t& dup);
  void RemoveDuplicates();
  void RemoveNameTitle();

protected:
  static constexpr double init_Ekin_nominal = 0;  //! Do not persistify!
  static constexpr double init_Ekin_min = 0;  //! Do not persistify!
  static constexpr double init_Ekin_max = 14000000;  //! Do not persistify!
  static constexpr double init_eta_nominal = 0;  //! Do not persistify!
  static constexpr double init_eta_min = -100;  //! Do not persistify!
  static constexpr double init_eta_max = 100;  //! Do not persistify!

  auto compare(const TFCSParametrizationBase& ref) const -> bool;

public:
  /// Update outputlevel
  /// for multiple levels
  using MLogging::setLevel;
  virtual void setLevel(int level, bool recursive)
  {
    this->setLevel(level);
    if (recursive)
      for (unsigned int i = 0; i < size(); ++i)
        (*this)[i]->setLevel(level, recursive);
  }

private:
  ClassDefOverride(TFCSParametrizationBase, 5)  // TFCSParametrizationBase
};
