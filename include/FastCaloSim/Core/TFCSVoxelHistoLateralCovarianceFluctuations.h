/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#ifndef TFCSVoxelHistoLateralCovarianceFluctuations_h
#define TFCSVoxelHistoLateralCovarianceFluctuations_h

#include "FastCaloSim/Core/TFCSLateralShapeParametrizationHitBase.h"
#include "FastCaloSim/Geometry/FastCaloSim_CaloCell_ID.h"
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
  virtual ~TFCSVoxelHistoLateralCovarianceFluctuations();

  virtual void set_geometry(CaloGeo* geo) override { m_geo = geo; };
  CaloGeo* get_geometry() { return m_geo; };

  bool initialize(TFile* inputfile, const std::string& folder);

  /// create one fluctuated shape for a shower to be applied as scale factor to
  /// the average shape Store the fluctuation in simulstate
  virtual FCSReturnCode simulate(
      TFCSSimulationState& simulstate,
      const TFCSTruthState* truth,
      const TFCSExtrapolationState* extrapol) const override;

  /// weight the energy of one hit by the fluctuation calculated in
  /// simulate(...)
  virtual FCSReturnCode simulate_hit(
      Hit& hit,
      TFCSSimulationState& simulstate,
      const TFCSTruthState* truth,
      const TFCSExtrapolationState* extrapol) override;

protected:
  CaloGeo* m_geo;  //! do not persistify

  static const std::uint32_t
      s_layer_hash[CaloCell_ID_FCS::MaxSample];  //! do not persistify
  static const std::uint32_t
      s_layer_hash_geo[CaloCell_ID_FCS::MaxSample];  //! do not persistify

  typedef std::vector<std::vector<float>> weight_t;

  std::vector<std::vector<std::vector<TFCS1DFunction*>>> m_transform;

  // For a 5*5 cell grid, nDim should be 5
  int m_nDim_x;
  int m_nDim_y;
  std::vector<TH2*> m_voxel_template;
  std::vector<TVectorD> m_parMeans;
  std::vector<TMatrixD> m_EigenVectors;  // Eigen-vectors of covariance
  std::vector<TVectorD> m_EigenValues;  // Eigen-values of covariance

  void MultiGaus(TFCSSimulationState& simulstate, TVectorD& genPars) const;

  ClassDefOverride(TFCSVoxelHistoLateralCovarianceFluctuations,
                   1)  // TFCSVoxelHistoLateralCovarianceFluctuations
};

#endif
