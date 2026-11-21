// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#pragma once

#include <map>
#include <mutex>
#include <set>
#include <string>

#include <FastCaloSim/FastCaloSim_export.h>

#include "FastCaloSim/Core/FCSReturnCode.h"
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

class FASTCALOSIM_EXPORT TFCSParametrizationBase
    : public TNamed
    , public FastCaloSim::MLogging
{
public:
  explicit TFCSParametrizationBase(const char* name = nullptr,
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
  }
  virtual auto is_match_Ekin(float /*Ekin*/) const -> bool { return false; }
  virtual auto is_match_eta(float /*eta*/) const -> bool { return false; }

  virtual auto is_match_Ekin_bin(std::size_t /*Ekin_bin*/) const -> bool
  {
    return false;
  }
  virtual auto is_match_calosample(std::size_t /*calosample*/) const -> bool
  {
    return false;
  }

  virtual auto is_match_all_pdgid() const -> bool
  {
    return TestBit(kMatchAllPDGID);
  }
  virtual auto is_match_all_Ekin() const -> bool { return false; }
  virtual auto is_match_all_eta() const -> bool { return false; }
  virtual auto is_match_all_Ekin_bin() const -> bool { return false; }
  virtual auto is_match_all_calosample() const -> bool { return false; }

  virtual auto pdgid() const -> const std::set<int>&
  {
    static const std::set<int> empty;
    return empty;
  }
  virtual auto Ekin_nominal() const -> double { return init_Ekin_nominal; }
  virtual auto Ekin_min() const -> double { return init_Ekin_min; }
  virtual auto Ekin_max() const -> double { return init_Ekin_max; }
  virtual auto eta_nominal() const -> double { return init_eta_nominal; }
  virtual auto eta_min() const -> double { return init_eta_min; }
  virtual auto eta_max() const -> double { return init_eta_max; }

  virtual void set_match_all_pdgid() { SetBit(kMatchAllPDGID); }
  virtual void reset_match_all_pdgid() { ResetBit(kMatchAllPDGID); }

  /// Method to set the geometry access pointer. Loops over daughter objects if
  /// present
  virtual void set_geometry(CaloGeo* geo);

  /// Some derived classes have daughter instances of TFCSParametrizationBase
  /// objects The size() and operator[] methods give general access to these
  /// daughters
  virtual auto size() const -> std::size_t { return 0; }

  /// Some derived classes have daughter instances of TFCSParametrizationBase
  /// objects The size() and operator[] methods give general access to these
  /// daughters
  virtual auto operator[](std::size_t /*ind*/) const
      -> const TFCSParametrizationBase*
  {
    return nullptr;
  }

  /// Some derived classes have daughter instances of TFCSParametrizationBase
  /// objects The size() and operator[] methods give general access to these
  /// daughters
  virtual auto operator[](std::size_t /*ind*/) -> TFCSParametrizationBase*
  {
    return nullptr;
  }

  /// Some derived classes have daughter instances of TFCSParametrizationBase
  /// objects The set_daughter method allows to change these daughters - expert
  /// use only! The original element at this position is not deleted
  virtual void set_daughter(std::size_t /*ind*/,
                            TFCSParametrizationBase* /*param*/)
  {
  }

  /// The == operator compares the content of instances.
  /// The implementation in the base class only returns true for a comparison
  /// with itself
  virtual auto operator==(const TFCSParametrizationBase& ref) const -> bool
  {
    return compare(ref);
  }

  /// Method in all derived classes to do some simulation
  virtual auto simulate(TFCSSimulationState& simulstate,
                        const TFCSTruthState* truth,
                        const TFCSExtrapolationState* extrapol) const
      -> FCSReturnCode;

  /// Called immediately after being read from file to fix histograms that
  /// erroneously have kIsOnHeap set. See ATLASSIM-7031.
  virtual void fixHists() {}

  /// Print object information.
  void Print(Option_t* option = "") const override;

  struct Duplicate_t
  {
    TFCSParametrizationBase* replace = nullptr;
    std::vector<TFCSParametrizationBase*> mother;
    std::vector<std::size_t> index;
  };
  using FindDuplicates_t = std::map<TFCSParametrizationBase*, Duplicate_t>;
  using FindDuplicateClasses_t = std::map<std::string, FindDuplicates_t>;
  void FindDuplicates(FindDuplicateClasses_t& dupclasses);
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

  /// Options of Print(): "short" anywhere in the option asks for a one-line
  /// summary, the rest of the option is the indentation prefix of every line.
  struct PrintOptions
  {
    bool shortprint;  ///< the option contains "short"
    bool longprint;  ///< print details: DEBUG level, or INFO without "short"
    std::string optprint;  ///< the option with every "short" removed
  };
  auto print_options(Option_t* option) const -> PrintOptions;

public:
  /// Update outputlevel
  /// for multiple levels
  using MLogging::setLevel;
  virtual void setLevel(int level, bool recursive)
  {
    this->setLevel(level);
    if (recursive) {
      for (std::size_t i = 0; i < size(); ++i) {
        if (TFCSParametrizationBase* param = (*this)[i]) {
          param->setLevel(level, recursive);
        }
      }
    }
  }

private:
  ClassDefOverride(TFCSParametrizationBase, 5)  // TFCSParametrizationBase
};
