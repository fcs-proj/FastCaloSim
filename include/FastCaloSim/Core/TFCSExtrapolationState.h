// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#pragma once

#include <cstddef>
#include <map>
#include <utility>

#include <FastCaloSim/FastCaloSim_export.h>
#include <TObject.h>

#include "FastCaloSim/Core/MLogging.h"
#include "FastCaloSim/Geometry/SubPos.h"

class FASTCALOSIM_EXPORT TFCSExtrapolationState
    : public TObject
    , public FastCaloSim::MLogging
{
public:
  TFCSExtrapolationState();

  void clear();

  static constexpr int NumSubPos = 3;

  void set_OK(std::size_t layer, FastCaloSim::SubPos subpos, bool val = true)
  {
    m_CaloOK[{layer, subpos}] = val;
  }

  void set_eta(std::size_t layer, FastCaloSim::SubPos subpos, double val)
  {
    m_etaCalo[{layer, subpos}] = val;
  }

  void set_phi(std::size_t layer, FastCaloSim::SubPos subpos, double val)
  {
    m_phiCalo[{layer, subpos}] = val;
  }

  void set_r(std::size_t layer, FastCaloSim::SubPos subpos, double val)
  {
    m_rCalo[{layer, subpos}] = val;
  }

  void set_z(std::size_t layer, FastCaloSim::SubPos subpos, double val)
  {
    m_zCalo[{layer, subpos}] = val;
  }

  void set_d(std::size_t layer, FastCaloSim::SubPos subpos, double val)
  {
    m_dCalo[{layer, subpos}] = val;
  }

  void set_detaBorder(std::size_t layer, FastCaloSim::SubPos subpos, double val)
  {
    m_distetaCaloBorder[{layer, subpos}] = val;
  }

  void set_IDCaloBoundary_eta(double val) { m_IDCaloBoundary_eta = val; }
  void set_IDCaloBoundary_phi(double val) { m_IDCaloBoundary_phi = val; }
  void set_IDCaloBoundary_r(double val) { m_IDCaloBoundary_r = val; }
  void set_IDCaloBoundary_x(double val) { m_IDCaloBoundary_x = val; }
  void set_IDCaloBoundary_y(double val) { m_IDCaloBoundary_y = val; }
  void set_IDCaloBoundary_z(double val) { m_IDCaloBoundary_z = val; }

  auto OK(std::size_t layer, FastCaloSim::SubPos subpos) const -> bool
  {
    // layers/subpositions never reached by the extrapolation are not stored;
    // read as the pre-map default
    auto it = m_CaloOK.find({layer, subpos});
    return it != m_CaloOK.end() ? it->second : false;
  }

  auto eta(std::size_t layer, FastCaloSim::SubPos subpos) const -> double
  {
    // layers/subpositions never reached by the extrapolation are not stored;
    // read as the pre-map default
    auto it = m_etaCalo.find({layer, subpos});
    return it != m_etaCalo.end() ? it->second : -999;
  }

  auto phi(std::size_t layer, FastCaloSim::SubPos subpos) const -> double
  {
    // layers/subpositions never reached by the extrapolation are not stored;
    // read as the pre-map default
    auto it = m_phiCalo.find({layer, subpos});
    return it != m_phiCalo.end() ? it->second : -999;
  }

  auto r(std::size_t layer, FastCaloSim::SubPos subpos) const -> double
  {
    // layers/subpositions never reached by the extrapolation are not stored;
    // read as the pre-map default
    auto it = m_rCalo.find({layer, subpos});
    return it != m_rCalo.end() ? it->second : 0;
  }

  auto z(std::size_t layer, FastCaloSim::SubPos subpos) const -> double
  {
    // layers/subpositions never reached by the extrapolation are not stored;
    // read as the pre-map default
    auto it = m_zCalo.find({layer, subpos});
    return it != m_zCalo.end() ? it->second : 0;
  }

  auto d(std::size_t layer, FastCaloSim::SubPos subpos) const -> double
  {
    // layers/subpositions never reached by the extrapolation are not stored;
    // read as the pre-map default
    auto it = m_dCalo.find({layer, subpos});
    return it != m_dCalo.end() ? it->second : 0;
  }

  auto detaBorder(std::size_t layer, FastCaloSim::SubPos subpos) const -> double
  {
    // layers/subpositions never reached by the extrapolation are not stored;
    // read as the pre-map default
    auto it = m_distetaCaloBorder.find({layer, subpos});
    return it != m_distetaCaloBorder.end() ? it->second : 0;
  }

  // Boundary getters
  auto IDCaloBoundary_eta() const -> double { return m_IDCaloBoundary_eta; }
  auto IDCaloBoundary_phi() const -> double { return m_IDCaloBoundary_phi; }
  auto IDCaloBoundary_r() const -> double { return m_IDCaloBoundary_r; }
  auto IDCaloBoundary_x() const -> double { return m_IDCaloBoundary_x; }
  auto IDCaloBoundary_y() const -> double { return m_IDCaloBoundary_y; }
  auto IDCaloBoundary_z() const -> double { return m_IDCaloBoundary_z; }

  auto IDCaloBoundary_AngleEta() const -> double
  {
    return m_IDCaloBoundary_AngleEta;
  }
  auto IDCaloBoundary_Angle3D() const -> double
  {
    return m_IDCaloBoundary_Angle3D;
  }

  // Simplified setting for angles
  void set_IDCaloBoundary_AngleEta(double val)
  {
    m_IDCaloBoundary_AngleEta = val;
  }

  void set_IDCaloBoundary_Angle3D(double val)
  {
    m_IDCaloBoundary_Angle3D = val;
  }
  void Print(Option_t* option = "") const override;

private:
  // 2D Map: (layer, subpos) -> extrapolated position
  using Key = std::pair<std::size_t, FastCaloSim::SubPos>;
  std::map<Key, bool> m_CaloOK;
  std::map<Key, double> m_etaCalo;
  std::map<Key, double> m_phiCalo;
  std::map<Key, double> m_rCalo;
  std::map<Key, double> m_zCalo;
  std::map<Key, double> m_dCalo;
  std::map<Key, double> m_distetaCaloBorder;

  double m_IDCaloBoundary_eta {-999};
  double m_IDCaloBoundary_phi {-999};
  double m_IDCaloBoundary_r {0};
  double m_IDCaloBoundary_x {0};
  double m_IDCaloBoundary_y {0};
  double m_IDCaloBoundary_z {0};

  double m_IDCaloBoundary_AngleEta {-999};
  double m_IDCaloBoundary_Angle3D {-999};

  ClassDefOverride(TFCSExtrapolationState, 3)  // TFCSExtrapolationState
};
