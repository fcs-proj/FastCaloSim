/*
  Copyright (C) 2002-2023 CERN for the benefit of the ATLAS collaboration
*/

#ifndef ISF_FASTCALOSIMEVENT_TFCSParametrizationBase_h
#define ISF_FASTCALOSIMEVENT_TFCSParametrizationBase_h

#include <FastCaloSim/FastCaloSim_export.h>

#include <map>
#include <mutex>
#include <set>

#include "FastCaloSim/MLogging.h"

class ICaloGeometry;
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
    , public ISF_FCS::MLogging
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

  virtual bool is_match_pdgid(int /*id*/) const
  {
    return TestBit(kMatchAllPDGID);
  };
  virtual bool is_match_Ekin(float /*Ekin*/) const { return false; };
  virtual bool is_match_eta(float /*eta*/) const { return false; };

  virtual bool is_match_Ekin_bin(int /*Ekin_bin*/) const { return false; };
  virtual bool is_match_calosample(int /*calosample*/) const { return false; };

  virtual bool is_match_all_pdgid() const { return TestBit(kMatchAllPDGID); };
  virtual bool is_match_all_Ekin() const { return false; };
  virtual bool is_match_all_eta() const { return false; };
  virtual bool is_match_all_Ekin_bin() const { return false; };
  virtual bool is_match_all_calosample() const { return false; };

  virtual const std::set<int>& pdgid() const
  {
    static const std::set<int> empty;
    return empty;
  };
  virtual double Ekin_nominal() const { return init_Ekin_nominal; };
  virtual double Ekin_min() const { return init_Ekin_min; };
  virtual double Ekin_max() const { return init_Ekin_max; };
  virtual double eta_nominal() const { return init_eta_nominal; };
  virtual double eta_min() const { return init_eta_min; };
  virtual double eta_max() const { return init_eta_max; };

  virtual void set_match_all_pdgid() { SetBit(kMatchAllPDGID); };
  virtual void reset_match_all_pdgid() { ResetBit(kMatchAllPDGID); };

  /// Method to set the geometry access pointer. Loops over daughter objects if
  /// present
  virtual void set_geometry(ICaloGeometry* geo);

  /// Some derived classes have daughter instances of TFCSParametrizationBase
  /// objects The size() and operator[] methods give general access to these
  /// daughters
  virtual unsigned int size() const { return 0; };

  /// Some derived classes have daughter instances of TFCSParametrizationBase
  /// objects The size() and operator[] methods give general access to these
  /// daughters
  virtual const TFCSParametrizationBase* operator[](unsigned int /*ind*/) const
  {
    return nullptr;
  };

  /// Some derived classes have daughter instances of TFCSParametrizationBase
  /// objects The size() and operator[] methods give general access to these
  /// daughters
  virtual TFCSParametrizationBase* operator[](unsigned int /*ind*/)
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
  virtual bool operator==(const TFCSParametrizationBase& ref) const
  {
    return compare(ref);
  };

  /// Method in all derived classes to do some simulation
  virtual FCSReturnCode simulate(TFCSSimulationState& simulstate,
                                 const TFCSTruthState* truth,
                                 const TFCSExtrapolationState* extrapol) const;

  /// Method in all derived classes to delete objects stored in the simulstate
  /// AuxInfo
  virtual void CleanAuxInfo(TFCSSimulationState& /*simulstate*/) const {};

  /// Print object information.
  void Print(Option_t* option = "") const;

  struct Duplicate_t
  {
    TFCSParametrizationBase* replace = nullptr;
    std::vector<TFCSParametrizationBase*> mother;
    std::vector<unsigned int> index;
  };
  typedef std::map<TFCSParametrizationBase*, Duplicate_t> FindDuplicates_t;
  typedef std::map<std::string, FindDuplicates_t> FindDuplicateClasses_t;
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

  bool compare(const TFCSParametrizationBase& ref) const;

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
  ClassDef(TFCSParametrizationBase, 3)  // TFCSParametrizationBase
};

#endif  // End header guards
