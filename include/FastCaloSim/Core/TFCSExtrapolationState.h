// Copyright (c) 2024 CERN for the benefit of the FastCaloSim project

#ifndef ISF_FASTCALOSIMEVENT_TFCSExtrapolationState_h
#define ISF_FASTCALOSIMEVENT_TFCSExtrapolationState_h

#include <map>
#include <utility>

#include <FastCaloSim/FastCaloSim_export.h>
#include <TObject.h>

#include "FastCaloSim/Core/MLogging.h"

class FASTCALOSIM_EXPORT TFCSExtrapolationState
    : public TObject
    , public ISF_FCS::MLogging
{
public:
  TFCSExtrapolationState();

  void clear();

  static constexpr int NumSubPos = 3;

  void set_OK(int layer, int subpos, bool val = true)
  {
    m_CaloOK[{layer, subpos}] = val;
  }

  void set_eta(int layer, int subpos, double val)
  {
    m_etaCalo[{layer, subpos}] = val;
  }

  void set_phi(int layer, int subpos, double val)
  {
    m_phiCalo[{layer, subpos}] = val;
  }

  void set_r(int layer, int subpos, double val)
  {
    m_rCalo[{layer, subpos}] = val;
  }

  void set_z(int layer, int subpos, double val)
  {
    m_zCalo[{layer, subpos}] = val;
  }

  void set_d(int layer, int subpos, double val)
  {
    m_dCalo[{layer, subpos}] = val;
  }

  void set_detaBorder(int layer, int subpos, double val)
  {
    m_distetaCaloBorder[{layer, subpos}] = val;
  }

  void set_IDCaloBoundary_eta(double val) { m_IDCaloBoundary_eta = val; }
  void set_IDCaloBoundary_phi(double val) { m_IDCaloBoundary_phi = val; }
  void set_IDCaloBoundary_r(double val) { m_IDCaloBoundary_r = val; }
  void set_IDCaloBoundary_x(double val) { m_IDCaloBoundary_x = val; }
  void set_IDCaloBoundary_y(double val) { m_IDCaloBoundary_y = val; }
  void set_IDCaloBoundary_z(double val) { m_IDCaloBoundary_z = val; }

  auto OK(int layer, int subpos) const -> bool
  {
    // layers/subpositions never reached by the extrapolation are not stored;
    // read as the pre-map default
    auto it = m_CaloOK.find({layer, subpos});
    return it != m_CaloOK.end() ? it->second : false;
  }

  auto eta(int layer, int subpos) const -> double
  {
    // layers/subpositions never reached by the extrapolation are not stored;
    // read as the pre-map default
    auto it = m_etaCalo.find({layer, subpos});
    return it != m_etaCalo.end() ? it->second : -999;
  }

  auto phi(int layer, int subpos) const -> double
  {
    // layers/subpositions never reached by the extrapolation are not stored;
    // read as the pre-map default
    auto it = m_phiCalo.find({layer, subpos});
    return it != m_phiCalo.end() ? it->second : -999;
  }

  auto r(int layer, int subpos) const -> double
  {
    // layers/subpositions never reached by the extrapolation are not stored;
    // read as the pre-map default
    auto it = m_rCalo.find({layer, subpos});
    return it != m_rCalo.end() ? it->second : 0;
  }

  auto z(int layer, int subpos) const -> double
  {
    // layers/subpositions never reached by the extrapolation are not stored;
    // read as the pre-map default
    auto it = m_zCalo.find({layer, subpos});
    return it != m_zCalo.end() ? it->second : 0;
  }

  auto d(int layer, int subpos) const -> double
  {
    // layers/subpositions never reached by the extrapolation are not stored;
    // read as the pre-map default
    auto it = m_dCalo.find({layer, subpos});
    return it != m_dCalo.end() ? it->second : 0;
  }

  auto detaBorder(int layer, int subpos) const -> double
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
  void Print(Option_t* option = "") const;

private:
  // 2D Map: (layer, subpos) -> extrapolated position
  std::map<std::pair<int, int>, bool> m_CaloOK;
  std::map<std::pair<int, int>, double> m_etaCalo;
  std::map<std::pair<int, int>, double> m_phiCalo;
  std::map<std::pair<int, int>, double> m_rCalo;
  std::map<std::pair<int, int>, double> m_zCalo;
  std::map<std::pair<int, int>, double> m_dCalo;
  std::map<std::pair<int, int>, double> m_distetaCaloBorder;

  double m_IDCaloBoundary_eta;
  double m_IDCaloBoundary_phi;
  double m_IDCaloBoundary_r;
  double m_IDCaloBoundary_x;
  double m_IDCaloBoundary_y;
  double m_IDCaloBoundary_z;

  double m_IDCaloBoundary_AngleEta;
  double m_IDCaloBoundary_Angle3D;

  ClassDef(TFCSExtrapolationState, 2)  // TFCSExtrapolationState
};

#endif
