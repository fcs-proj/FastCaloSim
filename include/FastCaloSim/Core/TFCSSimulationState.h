// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#ifndef ISF_FASTCALOSIMEVENT_TFCSSimulationState_h
#define ISF_FASTCALOSIMEVENT_TFCSSimulationState_h

#include <cstdint>
#include <set>
#include <unordered_map>

#include <FastCaloSim/FastCaloSim_export.h>
#include <TObject.h>

#include "FastCaloSim/Core/MLogging.h"
class TFCSParametrizationBase;

namespace CLHEP
{
class HepRandomEngine;
}

constexpr auto operator""_FCShash(char const* s, std::size_t count)
    -> std::uint32_t;

class FASTCALOSIM_EXPORT TFCSSimulationState
    : public TObject
    , public ISF_FCS::MLogging
{
public:
  TFCSSimulationState(CLHEP::HepRandomEngine* randomEngine = nullptr);

  auto randomEngine() -> CLHEP::HepRandomEngine* { return m_randomEngine; }
  void setRandomEngine(CLHEP::HepRandomEngine* engine)
  {
    m_randomEngine = engine;
  }

  auto is_valid() const -> bool { return m_Ebin >= 0; };
  auto E() const -> double { return m_Etot; };
  // NOTE: layers without energy are not stored in the maps; as for the
  // original array-based implementation, such layers read as 0 energy.
  auto E(int sample) const -> double
  {
    auto it = m_E.find(sample);
    return it != m_E.end() ? it->second : 0.;
  };
  auto Efrac(int sample) const -> double
  {
    auto it = m_Efrac.find(sample);
    return it != m_Efrac.end() ? it->second : 0.;
  };
  auto Ebin() const -> int { return m_Ebin; };

  void set_Ebin(int bin) { m_Ebin = bin; };
  void set_E(int sample, double Esample) { m_E[sample] = Esample; };
  void set_Efrac(int sample, double Efracsample)
  {
    m_Efrac[sample] = Efracsample;
  };
  void set_E(double E) { m_Etot = E; };
  void add_E(int sample, double Esample)
  {
    m_E[sample] += Esample;
    m_Etot += Esample;
  };

  // maps the cell id to the energy deposited in the cell
  using cellmap = std::unordered_map<unsigned long long, float>;

  auto cells() -> cellmap& { return m_cells; };
  auto cells() const -> const cellmap& { return m_cells; };

  void deposit(const unsigned long long cell_id, float E);

  void Print(Option_t* option = "") const override;

  // TODO: Remove explicit functions for SF here and use
  // getAuxInfo<double>("SF"_FCShash) and setAuxInfo<double>("SF"_FCShash,mysf)
  // directly in the energy parametrization
  void set_SF(double mysf) { setAuxInfo<double>("SF"_FCShash, mysf); };
  auto get_SF() -> double { return getAuxInfo<double>("SF"_FCShash); }

  void clear();

private:
  CLHEP::HepRandomEngine* m_randomEngine;

  int m_Ebin;
  double m_Etot;
  // TO BE CLEANED UP! SHOULD ONLY STORE EITHER E OR EFRAC!!!
  std::unordered_map<int, double> m_E;
  std::unordered_map<int, double> m_Efrac;

  cellmap m_cells;

public:
  // Allow to store arbitrary type objects as auxiliary information
  // Use compile time hashes of strings as index to an unordered map of union
  // AuxInfo_t Example: TFCSSimulationState s;
  // s.setAuxInfo<double>("SF"_FCShash,2);
  //
  // If pointers are stored, a dedicated cleanup is needed
  // If a new data type is needed, a cast operator and an explicit template
  // implementation of the set method has to be added
  union AuxInfo_t
  {
    bool b;
    char c;
    int i;
    float f;
    double d;
    void* p;

    // cast operators
    operator bool() const { return b; };
    operator char() const { return c; };
    operator int() const { return i; };
    operator float() const { return f; };
    operator double() const { return d; };
    operator void*() const { return p; };

    // template set method. No general implementation exist, only explicit
    // implementations are added after the class definition
    template<class T>
    void set(T val);
    void set(const AuxInfo_t& val);
  };

  // FNV-1a 32bit hashing algorithm that is evaluated during compile time
  // function taken from https://gist.github.com/Lee-R/3839813
  static constexpr auto fnv1a_32(char const* s, std::size_t count)
      -> std::uint32_t
  {
    return ((count ? fnv1a_32(s, count - 1) : 2166136261u) ^ s[count])
        * 16777619u;
  }
  // Run time call for hash function
  static auto getAuxIndex(const std::string& s) -> std::uint32_t;
  static auto getAuxIndex(const char* s) -> std::uint32_t;

  // Check if some auxiliary information is stored
  auto hasAuxInfo(std::uint32_t index) const -> bool
  {
    return m_AuxInfo.find(index) != m_AuxInfo.end();
  };

  // Get auxiliary info
  // Use as TFCSSimulationState::getAuxInfo<int>(index)
  template<class T>
  inline auto getAuxInfo(std::uint32_t index) const -> const T
  {
    return static_cast<T>(m_AuxInfo.at(index));
  }

  // Set auxiliary info
  // Use as TFCSSimulationState::setAuxInfo<double>(7,2.0f)
  //     or TFCSSimulationState::setAuxInfo(7,2.0)
  template<class T>
  inline void setAuxInfo(std::uint32_t index, const T& val)
  {
    m_AuxInfo[index].set<T>(val);
  }

  void AddAuxInfoCleanup(const TFCSParametrizationBase* para);
  void DoAuxInfoCleanup();

private:
  std::unordered_map<std::uint32_t, AuxInfo_t> m_AuxInfo;  //! Do not persistify
  std::set<const TFCSParametrizationBase*>
      m_AuxInfoCleanup;  //! Do not persistify

  ClassDef(TFCSSimulationState, 3)  // TFCSSimulationState
};

// Explicit template implementations for template<class T> void
// TFCSSimulationState::AuxInfo_t::set(T val);
inline void TFCSSimulationState::AuxInfo_t::set(
    const TFCSSimulationState::AuxInfo_t& val)
{
  *this = val;
}
template<>
inline void TFCSSimulationState::AuxInfo_t::set<bool>(bool val)
{
  b = val;
}
template<>
inline void TFCSSimulationState::AuxInfo_t::set<char>(char val)
{
  c = val;
}
template<>
inline void TFCSSimulationState::AuxInfo_t::set<int>(int val)
{
  i = val;
}
template<>
inline void TFCSSimulationState::AuxInfo_t::set<float>(float val)
{
  f = val;
}
template<>
inline void TFCSSimulationState::AuxInfo_t::set<double>(double val)
{
  d = val;
}
template<>
inline void TFCSSimulationState::AuxInfo_t::set<void*>(void* val)
{
  p = val;
}

// Implementation of the compile time text hash operator that can be used for
// human readable indices to the AuxInfo
constexpr auto operator""_FCShash(char const* s, std::size_t count)
    -> std::uint32_t
{
  return TFCSSimulationState::fnv1a_32(s, count);
}

#endif
