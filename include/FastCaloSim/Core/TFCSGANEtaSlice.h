// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

//////////////////////////////////////////////////////////////////
// TFCSGANEtaSlice.h, (c) ATLAS Detector software
///////////////////////////////////////////////////////////////////

#pragma once

// STL includes
#include <map>
#include <vector>

#include "FastCaloSim/Core/MLogging.h"
#include "FastCaloSim/Core/TFCSExtrapolationState.h"
#include "FastCaloSim/Core/TFCSGANXMLParameters.h"
#include "FastCaloSim/Core/TFCSSimulationState.h"
#include "FastCaloSim/Core/TFCSTruthState.h"

// generic network class
#include "FastCaloSim/Core/VNetworkBase.h"
// net class for legacy loading
#include <fstream>

#include "FastCaloSim/Core/TFCSGANLWTNNHandler.h"

class TFCSGANEtaSlice : public FastCaloSim::MLogging
{
public:
  TFCSGANEtaSlice();
  TFCSGANEtaSlice(int pid,
                  int etaMin,
                  int etaMax,
                  const TFCSGANXMLParameters& param);
  ~TFCSGANEtaSlice() override;

  TFCSGANEtaSlice(const TFCSGANEtaSlice&) = delete;
  auto operator=(const TFCSGANEtaSlice&) -> TFCSGANEtaSlice& = delete;

  using FitResultsPerLayer = std::map<int, std::vector<double>>;
  using ExtrapolatorWeights = std::map<int, double>;
  using NetworkInputs = std::map<std::string, std::map<std::string, double>>;
  using NetworkOutputs = std::map<std::string, double>;

  auto LoadGAN() -> bool;
  void CalculateMeanPointFromDistributionOfR();
  void ExtractExtrapolatorMeansFromInputs();

  auto GetNetworkOutputs(const TFCSTruthState* truth,
                         const TFCSExtrapolationState* extrapol,
                         TFCSSimulationState& simulstate) const
      -> NetworkOutputs;

  auto IsGanCorrectlyLoaded() const -> bool;
  auto GetFitResults() const -> const FitResultsPerLayer&
  {
    return m_allFitResults;
  }
  auto GetExtrapolatorWeights() const -> const ExtrapolatorWeights&
  {
    return m_extrapolatorWeights;
  }

  void Print() const;

private:
  int m_pid {};
  int m_etaMin {};
  int m_etaMax {};

  std::string m_inputFolderName;

  FitResultsPerLayer m_allFitResults;
  ExtrapolatorWeights m_extrapolatorWeights;

  // legacy - keep or streamers are confused by
  // old classes that didn't inherit
  TFCSGANLWTNNHandler* m_gan_all = nullptr;
  TFCSGANLWTNNHandler* m_gan_low = nullptr;
  TFCSGANLWTNNHandler* m_gan_high = nullptr;
  // updated - can take an old or new class
  std::unique_ptr<VNetworkBase> m_net_all = nullptr;
  std::unique_ptr<VNetworkBase> m_net_low = nullptr;
  std::unique_ptr<VNetworkBase> m_net_high = nullptr;
  // getters so that we are insensitive to where the data actually is
  auto GetNetAll() const -> VNetworkBase*;
  auto GetNetLow() const -> VNetworkBase*;
  auto GetNetHigh() const -> VNetworkBase*;

  auto LoadGANNoRange(std::string inputFileName) -> bool;
  auto LoadGANFromRange(std::string inputFileName, std::string energyRange)
      -> bool;

  TFCSGANXMLParameters m_param;

  ClassDefOverride(TFCSGANEtaSlice, 6)  // TFCSGANEtaSlice
};
