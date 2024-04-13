/*
  Copyright (C) 2002-2023 CERN for the benefit of the ATLAS collaboration
*/

#ifndef ISF_FASTCALOSIMEVENT_TFCSSimulationState_h
#define ISF_FASTCALOSIMEVENT_TFCSSimulationState_h

#include <TObject.h>
#include "FastCaloSim/FastCaloSim_CaloCell_ID.h"
#include "FastCaloSim/MLogging.h"
#include <map>
#include <set>
#include <unordered_map>
#include <vector>
#include <cstdint>

#undef FCS_USE_HASH_SORTED_CELLMAP
#ifdef FCS_USE_HASH_SORTED_CELLMAP
#include "FastCaloSim/CaloDetDescrElement.h"
#else
class CaloDetDescrElement;
#endif
class TFCSParametrizationBase;

namespace CLHEP {
class HepRandomEngine;
}

constexpr std::uint32_t operator"" _FCShash(char const *s, std::size_t count);


class TFCSSimulationState : public TObject, public ISF_FCS::MLogging {
public:
  TFCSSimulationState(CLHEP::HepRandomEngine *randomEngine = nullptr);

  CLHEP::HepRandomEngine *randomEngine() { return m_randomEngine; }
  void setRandomEngine(CLHEP::HepRandomEngine *engine) {
    m_randomEngine = engine;
  }

  bool is_valid() const { return m_Ebin >= 0; };
  double E() const { return m_Etot; };
  double E(int sample) const { return m_E[sample]; };
  double Efrac(int sample) const { return m_Efrac[sample]; };
  int Ebin() const { return m_Ebin; };

  void set_Ebin(int bin) { m_Ebin = bin; };
  void set_E(int sample, double Esample) { m_E[sample] = Esample; };
  void set_Efrac(int sample, double Efracsample) {
    m_Efrac[sample] = Efracsample;
  };
  void set_E(double E) { m_Etot = E; };
  void add_E(int sample, double Esample) {
    m_E[sample] += Esample;
    m_Etot += Esample;
  };

#ifdef FCS_USE_HASH_SORTED_CELLMAP
  struct hashesCmp {
    bool operator()(const CaloDetDescrElement *a,
                    const CaloDetDescrElement *b) const {
      return a->calo_hash() < b->calo_hash();
    }
  };
  // Being able to force the order iteration over the Cellmap_t is very useful
  // when debugging small differences in output
  typedef std::map<const CaloDetDescrElement *, float, hashesCmp> Cellmap_t;
#else
  typedef std::map<const CaloDetDescrElement *, float> Cellmap_t;
#endif

  Cellmap_t &cells() { return m_cells; };
  const Cellmap_t &cells() const { return m_cells; };
  void deposit(const CaloDetDescrElement *cellele, float E);

  void Print(Option_t *option = "") const;

  // TODO: Remove explicit functions for SF here and use
  // getAuxInfo<double>("SF"_FCShash) and setAuxInfo<double>("SF"_FCShash,mysf)
  // directly in the energy parametrization
  void set_SF(double mysf) { setAuxInfo<double>("SF"_FCShash, mysf); };
  double get_SF() { return getAuxInfo<double>("SF"_FCShash); }

  void clear();

private:

  CLHEP::HepRandomEngine *m_randomEngine;

  int m_Ebin;
  double m_Etot;
  // TO BE CLEANED UP! SHOULD ONLY STORE EITHER E OR EFRAC!!!
  double m_E[CaloCell_ID_FCS::MaxSample];
  double m_Efrac[CaloCell_ID_FCS::MaxSample];

  Cellmap_t m_cells;

public:
  // Allow to store arbitrary type objects as auxilliary information
  // Use compile time hashes of strings as index to an unordered map of union
  // AuxInfo_t Example: TFCSSimulationState s;
  // s.setAuxInfo<double>("SF"_FCShash,2);
  //
  // If pointers are stored, a dedicated cleanup is needed
  // If a new data type is needed, a cast operator and an explicit template
  // implementation of the set method has to be added
  union AuxInfo_t {
    bool b;
    char c;
    int i;
    float f;
    double d;
    void *p;

    // cast operators
    operator bool() const { return b; };
    operator char() const { return c; };
    operator int() const { return i; };
    operator float() const { return f; };
    operator double() const { return d; };
    operator void *() const { return p; };

    // template set method. No general implementation exist, only explict
    // implementations are added after the class definition
    template <class T> void set(T val);
    void set(const AuxInfo_t &val);
  };

  // FNV-1a 32bit hashing algorithm that is evaluated during compile time
  // function taken from https://gist.github.com/Lee-R/3839813
  static constexpr std::uint32_t fnv1a_32(char const *s, std::size_t count) {
    return ((count ? fnv1a_32(s, count - 1) : 2166136261u) ^ s[count]) *
           16777619u;
  }
  // Run time call for hash function
  static std::uint32_t getAuxIndex(const std::string &s);
  static std::uint32_t getAuxIndex(const char *s);

  // Check if some auxiliary information is stored
  bool hasAuxInfo(std::uint32_t index) const {
    return m_AuxInfo.find(index) != m_AuxInfo.end();
  };

  // Get auxiliary info
  // Use as TFCSSimulationState::getAuxInfo<int>(index)
  template <class T> inline const T getAuxInfo(std::uint32_t index) const {
    return static_cast<T>(m_AuxInfo.at(index));
  }

  // Set auxiliary info
  // Use as TFCSSimulationState::setAuxInfo<double>(7,2.0f)
  //     or TFCSSimulationState::setAuxInfo(7,2.0)
  template <class T> inline void setAuxInfo(std::uint32_t index, const T &val) {
    m_AuxInfo[index].set<T>(val);
  }

  void AddAuxInfoCleanup(const TFCSParametrizationBase *para);
  void DoAuxInfoCleanup();

private:
  std::unordered_map<std::uint32_t, AuxInfo_t> m_AuxInfo; //! Do not persistify
  std::set<const TFCSParametrizationBase *>
      m_AuxInfoCleanup; //! Do not persistify

  ClassDef(TFCSSimulationState, 3) // TFCSSimulationState
};

// Explicit template implementations for template<class T> void
// TFCSSimulationState::AuxInfo_t::set(T val);
inline void
TFCSSimulationState::AuxInfo_t::set(const TFCSSimulationState::AuxInfo_t &val) {
  *this = val;
}
template <> inline void TFCSSimulationState::AuxInfo_t::set<bool>(bool val) {
  b = val;
}
template <> inline void TFCSSimulationState::AuxInfo_t::set<char>(char val) {
  c = val;
}
template <> inline void TFCSSimulationState::AuxInfo_t::set<int>(int val) {
  i = val;
}
template <> inline void TFCSSimulationState::AuxInfo_t::set<float>(float val) {
  f = val;
}
template <>
inline void TFCSSimulationState::AuxInfo_t::set<double>(double val) {
  d = val;
}
template <> inline void TFCSSimulationState::AuxInfo_t::set<void *>(void *val) {
  p = val;
}

// Implementation of the complile time text hash operator that can be used for
// human readable indices to the AuxInfo
constexpr std::uint32_t operator"" _FCShash(char const *s, std::size_t count) {
  return TFCSSimulationState::fnv1a_32(s, count);
}

#endif
