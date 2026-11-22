// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#pragma once

#include <any>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <unordered_map>
#include <utility>
#include <variant>

#include <FastCaloSim/FastCaloSim_export.h>
#include <TObject.h>

#include "FastCaloSim/Core/MLogging.h"

namespace CLHEP
{
class HepRandomEngine;
}  // namespace CLHEP

constexpr auto operator""_FCShash(char const* s, std::size_t count)
    -> std::uint32_t;

class FASTCALOSIM_EXPORT TFCSSimulationState
    : public TObject
    , public FastCaloSim::MLogging
{
public:
  explicit TFCSSimulationState(CLHEP::HepRandomEngine* randomEngine = nullptr);

  auto randomEngine() -> CLHEP::HepRandomEngine* { return m_randomEngine; }
  void setRandomEngine(CLHEP::HepRandomEngine* engine)
  {
    m_randomEngine = engine;
  }

  auto is_valid() const -> bool { return m_Ebin.has_value(); }
  auto E() const -> double { return m_Etot; }
  // NOTE: layers without energy are not stored in the maps; as for the
  // original array-based implementation, such layers read as 0 energy.
  auto E(std::size_t sample) const -> double
  {
    auto it = m_E.find(sample);
    return it != m_E.end() ? it->second : 0.;
  }
  auto Efrac(std::size_t sample) const -> double
  {
    auto it = m_Efrac.find(sample);
    return it != m_Efrac.end() ? it->second : 0.;
  }
  /// The energy bin of the PCA energy parametrization, starting at 1; no value
  /// before it is chosen
  auto Ebin() const -> std::optional<std::size_t> { return m_Ebin; }

  void set_Ebin(std::size_t bin) { m_Ebin = bin; }
  void set_E(std::size_t sample, double Esample) { m_E[sample] = Esample; }
  void set_Efrac(std::size_t sample, double Efracsample)
  {
    m_Efrac[sample] = Efracsample;
  }
  void set_E(double E) { m_Etot = E; }
  void add_E(std::size_t sample, double Esample)
  {
    m_E[sample] += Esample;
    m_Etot += Esample;
  }

  // maps the cell id to the energy deposited in the cell
  using cellmap = std::unordered_map<unsigned long long, double>;

  auto cells() -> cellmap& { return m_cells; }
  auto cells() const -> const cellmap& { return m_cells; }

  void deposit(unsigned long long cell_id, double E);

  void Print(Option_t* option = "") const override;

  // TODO: Remove explicit functions for SF here and use
  // getAuxInfo<double>("SF"_FCShash) and setAuxInfo<double>("SF"_FCShash,mysf)
  // directly in the energy parametrization
  void set_SF(double mysf) { setAuxInfo<double>("SF"_FCShash, mysf); }
  auto get_SF() const -> double { return getAuxInfo<double>("SF"_FCShash); }

  /// Reset the state for a new particle: energies, energy bin, deposited
  /// cells, auxiliary information and objects; the scale factor becomes 1.
  /// The random engine is kept.
  void clear();

private:
  CLHEP::HepRandomEngine* m_randomEngine;

  std::optional<std::size_t> m_Ebin;  //!
  double m_Etot {0};
  // TO BE CLEANED UP! SHOULD ONLY STORE EITHER E OR EFRAC!!!
  std::unordered_map<std::size_t, double> m_E;
  std::unordered_map<std::size_t, double> m_Efrac;

  cellmap m_cells;

public:
  // Auxiliary information: values that parametrizations store for each other
  // or for later steps of the same particle, indexed by compile time hashes
  // of names. Example: TFCSSimulationState s;
  // s.setAuxInfo<double>("SF"_FCShash, 2);
  //
  // Scalar values are stored with setAuxInfo<T>, objects (vectors, event
  // data, ...) with setAuxObject<T>, see below.
  /// The value of an auxiliary information entry: one of these types, the
  /// one it was stored with
  using AuxInfo_t = std::variant<bool, char, int, float, double>;

  // FNV-1a 32bit hashing algorithm that is evaluated during compile time
  // function taken from https://gist.github.com/Lee-R/3839813
  static constexpr auto fnv1a_32(char const* s, std::size_t count)
      -> std::uint32_t
  {
    return ((count != 0 ? fnv1a_32(s, count - 1) : 2166136261U)
            ^ static_cast<std::uint32_t>(s[count]))
        * 16777619U;
  }
  // Index of an auxiliary information entry per calorimeter layer: the hash of
  // a name continued with the layer number (one more FNV-1a step), so that
  // every name has its own indices for all layers.
  // Example: layerAuxIndex("ExtrapWeight"_FCShash, layer)
  static constexpr auto layerAuxIndex(std::uint32_t name_hash,
                                      std::size_t layer) -> std::uint32_t
  {
    return (name_hash ^ static_cast<std::uint32_t>(layer)) * 16777619U;
  }
  // Run time call for hash function
  static auto getAuxIndex(const std::string& s) -> std::uint32_t;
  static auto getAuxIndex(const char* s) -> std::uint32_t;

  // Check if some auxiliary information is stored
  auto hasAuxInfo(std::uint32_t index) const -> bool
  {
    return m_AuxInfo.contains(index);
  }

  // Get auxiliary info
  // Use as TFCSSimulationState::getAuxInfo<int>(index)
  template<class T>
  auto getAuxInfo(std::uint32_t index) const -> T
  {
    return std::get<T>(m_AuxInfo.at(index));
  }

  // Set auxiliary info
  // Use as TFCSSimulationState::setAuxInfo<double>(7,2.0f)
  //     or TFCSSimulationState::setAuxInfo(7,2.0)
  template<class T>
  void setAuxInfo(std::uint32_t index, const T& val)
  {
    m_AuxInfo[index] = val;
  }

  // Auxiliary objects: owned by the simulation state, so that they are
  // destroyed when they are replaced or erased, or with the state
  // Store an object, replacing a stored one, and return it
  template<class T>
  auto setAuxObject(std::uint32_t index, T object) -> T&
  {
    return m_AuxObjects[index].emplace<T>(std::move(object));
  }

  // The object stored at index, or nullptr if there is none of type T
  template<class T>
  auto getAuxObject(std::uint32_t index) -> T*
  {
    const auto it = m_AuxObjects.find(index);
    return it == m_AuxObjects.end() ? nullptr : std::any_cast<T>(&it->second);
  }
  template<class T>
  auto getAuxObject(std::uint32_t index) const -> const T*
  {
    const auto it = m_AuxObjects.find(index);
    return it == m_AuxObjects.end() ? nullptr : std::any_cast<T>(&it->second);
  }

  // Destroy the object stored at index, if any
  void eraseAuxObject(std::uint32_t index) { m_AuxObjects.erase(index); }

private:
  std::unordered_map<std::uint32_t, AuxInfo_t> m_AuxInfo;  //! Do not persistify
  std::unordered_map<std::uint32_t, std::any>
      m_AuxObjects;  //! Do not persistify

  ClassDefOverride(TFCSSimulationState, 5)  // TFCSSimulationState
};

// Implementation of the compile time text hash operator that can be used for
// human readable indices to the AuxInfo
constexpr auto operator""_FCShash(char const* s, std::size_t count)
    -> std::uint32_t
{
  return TFCSSimulationState::fnv1a_32(s, count);
}
