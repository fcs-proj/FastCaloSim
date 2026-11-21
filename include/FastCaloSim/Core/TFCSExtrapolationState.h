// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#pragma once

#include <cstddef>
#include <vector>

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

  static constexpr std::size_t NumSubPos = 3;

  void set_OK(std::size_t layer, FastCaloSim::SubPos subpos, bool val = true)
  {
    at(layer, subpos).ok = val;
  }

  void set_eta(std::size_t layer, FastCaloSim::SubPos subpos, double val)
  {
    at(layer, subpos).eta = val;
  }

  void set_phi(std::size_t layer, FastCaloSim::SubPos subpos, double val)
  {
    at(layer, subpos).phi = val;
  }

  void set_r(std::size_t layer, FastCaloSim::SubPos subpos, double val)
  {
    at(layer, subpos).r = val;
  }

  void set_z(std::size_t layer, FastCaloSim::SubPos subpos, double val)
  {
    at(layer, subpos).z = val;
  }

  void set_d(std::size_t layer, FastCaloSim::SubPos subpos, double val)
  {
    at(layer, subpos).d = val;
  }

  void set_detaBorder(std::size_t layer, FastCaloSim::SubPos subpos, double val)
  {
    at(layer, subpos).deta_border = val;
  }

  void set_IDCaloBoundary_eta(double val) { m_IDCaloBoundary_eta = val; }
  void set_IDCaloBoundary_phi(double val) { m_IDCaloBoundary_phi = val; }
  void set_IDCaloBoundary_r(double val) { m_IDCaloBoundary_r = val; }
  void set_IDCaloBoundary_x(double val) { m_IDCaloBoundary_x = val; }
  void set_IDCaloBoundary_y(double val) { m_IDCaloBoundary_y = val; }
  void set_IDCaloBoundary_z(double val) { m_IDCaloBoundary_z = val; }

  auto OK(std::size_t layer, FastCaloSim::SubPos subpos) const -> bool
  {
    return get(layer, subpos).ok;
  }

  auto eta(std::size_t layer, FastCaloSim::SubPos subpos) const -> double
  {
    return get(layer, subpos).eta;
  }

  auto phi(std::size_t layer, FastCaloSim::SubPos subpos) const -> double
  {
    return get(layer, subpos).phi;
  }

  auto r(std::size_t layer, FastCaloSim::SubPos subpos) const -> double
  {
    return get(layer, subpos).r;
  }

  auto z(std::size_t layer, FastCaloSim::SubPos subpos) const -> double
  {
    return get(layer, subpos).z;
  }

  auto d(std::size_t layer, FastCaloSim::SubPos subpos) const -> double
  {
    return get(layer, subpos).d;
  }

  auto detaBorder(std::size_t layer, FastCaloSim::SubPos subpos) const -> double
  {
    return get(layer, subpos).deta_border;
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

  /// Extrapolated position of one sub-position of a layer. The defaults are
  /// the values of a sub-position that the extrapolation did not reach.
  struct LayerPosition
  {
    bool ok {false};
    double eta {-999};
    double phi {-999};
    double r {0};
    double z {0};
    double d {0};
    double deta_border {0};
  };

private:
  static constexpr auto index(std::size_t layer, FastCaloSim::SubPos subpos)
      -> std::size_t
  {
    return (layer * NumSubPos) + static_cast<std::size_t>(subpos);
  }
  /// The stored position, added (with the defaults) if it is not stored yet
  auto at(std::size_t layer, FastCaloSim::SubPos subpos) -> LayerPosition&
  {
    const std::size_t i = index(layer, subpos);
    if (i >= m_positions.size()) {
      m_positions.resize((layer + 1) * NumSubPos);
    }
    return m_positions[i];
  }
  /// The stored position, or the defaults if it is not stored
  auto get(std::size_t layer, FastCaloSim::SubPos subpos) const
      -> const LayerPosition&
  {
    static constexpr LayerPosition not_reached {};
    const std::size_t i = index(layer, subpos);
    return i < m_positions.size() ? m_positions[i] : not_reached;
  }

  /// Positions of all layers, NumSubPos per layer (index(layer, subpos))
  std::vector<LayerPosition> m_positions;

  double m_IDCaloBoundary_eta {-999};
  double m_IDCaloBoundary_phi {-999};
  double m_IDCaloBoundary_r {0};
  double m_IDCaloBoundary_x {0};
  double m_IDCaloBoundary_y {0};
  double m_IDCaloBoundary_z {0};

  double m_IDCaloBoundary_AngleEta {-999};
  double m_IDCaloBoundary_Angle3D {-999};

  ClassDefOverride(TFCSExtrapolationState, 4)  // TFCSExtrapolationState
};
