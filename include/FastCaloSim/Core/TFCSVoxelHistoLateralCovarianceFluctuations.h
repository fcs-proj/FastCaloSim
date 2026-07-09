// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#pragma once

#include <string>

#include "FastCaloSim/Core/TFCSLateralShapeParametrizationHitBase.h"
#include "TH2.h"
#include "TMatrixD.h"
#include "TVectorD.h"

class TFCS1DFunction;
class TFile;

class TFCSVoxelHistoLateralCovarianceFluctuations
    : public TFCSLateralShapeParametrizationHitBase
{
public:
  TFCSVoxelHistoLateralCovarianceFluctuations(const char* name = nullptr,
                                              const char* title = nullptr);
  ~TFCSVoxelHistoLateralCovarianceFluctuations() override;

  void set_geometry(CaloGeo* geo) override { m_geo = geo; };
  auto get_geometry() -> CaloGeo* { return m_geo; };

  auto initialize(TFile* inputfile, const std::string& folder) -> bool;

  /// create one fluctuated shape for a shower to be applied as scale factor to
  /// the average shape Store the fluctuation in simulstate
  auto simulate(TFCSSimulationState& simulstate,
                const TFCSTruthState* truth,
                const TFCSExtrapolationState* extrapol) const
      -> FCSReturnCode override;

  /// weight the energy of one hit by the fluctuation calculated in
  /// simulate(...)
  auto simulate_hit(Hit& hit,
                    TFCSSimulationState& simulstate,
                    const TFCSTruthState* truth,
                    const TFCSExtrapolationState* extrapol)
      -> FCSReturnCode override;

protected:
  CaloGeo* m_geo = nullptr;  //! do not persistify

  using weight_t = std::vector<std::vector<float>>;

  std::vector<std::vector<std::vector<TFCS1DFunction*>>> m_transform;

  // For a 5*5 cell grid, nDim should be 5
  int m_nDim_x {};
  int m_nDim_y {};
  std::vector<TH2*> m_voxel_template;
  std::vector<TVectorD> m_parMeans;
  std::vector<TMatrixD> m_EigenVectors;  // Eigen-vectors of covariance
  std::vector<TVectorD> m_EigenValues;  // Eigen-values of covariance

  void MultiGaus(TFCSSimulationState& simulstate, TVectorD& genPars) const;

  ClassDefOverride(TFCSVoxelHistoLateralCovarianceFluctuations,
                   1)  // TFCSVoxelHistoLateralCovarianceFluctuations
};
