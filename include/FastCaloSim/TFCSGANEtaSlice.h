/*
  Copyright (C) 2002-2024 CERN for the benefit of the ATLAS collaboration
*/

//////////////////////////////////////////////////////////////////
// TFCSGANEtaSlice.h, (c) ATLAS Detector software
///////////////////////////////////////////////////////////////////

#ifndef ISF_TFCSGANETASLICE_H
#define ISF_TFCSGANETASLICE_H 1

// ISF includes
#include <map>
#include <vector>

#include "FastCaloSim/TFCSTruthState.h"
#include "FastCaloSim/TFCSSimulationState.h"
#include "FastCaloSim/TFCSExtrapolationState.h"
#include "FastCaloSim/TFCSGANXMLParameters.h"
#include "FastCaloSim/MLogging.h"

// generic network class
#include "FastCaloSim/VNetworkBase.h"
// net class for legacy loading
#include "FastCaloSim/TFCSGANLWTNNHandler.h"

#include <fstream>


class TFCSGANEtaSlice : public ISF_FCS::MLogging {
public:
  TFCSGANEtaSlice();
  TFCSGANEtaSlice(int pid, int etaMin, int etaMax,
                  const TFCSGANXMLParameters &param);
  virtual ~TFCSGANEtaSlice();

  typedef std::map<int, std::vector<double>> FitResultsPerLayer;
  typedef std::map<int, double> ExtrapolatorWeights;
  typedef std::map<std::string, std::map<std::string, double>> NetworkInputs;
  typedef std::map<std::string, double> NetworkOutputs;

  bool LoadGAN();
  void CalculateMeanPointFromDistributionOfR();
  void ExtractExtrapolatorMeansFromInputs();

  NetworkOutputs GetNetworkOutputs(const TFCSTruthState *truth,
                                   const TFCSExtrapolationState *extrapol,
                                   TFCSSimulationState simulstate) const;

  bool IsGanCorrectlyLoaded() const;
  FitResultsPerLayer GetFitResults() const { return m_allFitResults; }
  ExtrapolatorWeights GetExtrapolatorWeights() { return m_extrapolatorWeights; }

  void Print() const;

private:
  int m_pid;
  int m_etaMin;
  int m_etaMax;

  std::string m_inputFolderName;

  FitResultsPerLayer m_allFitResults;
  ExtrapolatorWeights m_extrapolatorWeights;

  // legacy - keep or streamers are confused by
  // old classes that didn't inherit
  TFCSGANLWTNNHandler *m_gan_all = nullptr;
  TFCSGANLWTNNHandler *m_gan_low = nullptr;
  TFCSGANLWTNNHandler *m_gan_high = nullptr;
  // updated - can take an old or new class
  std::unique_ptr<VNetworkBase> m_net_all = nullptr;
  std::unique_ptr<VNetworkBase> m_net_low = nullptr;
  std::unique_ptr<VNetworkBase> m_net_high = nullptr;
  // getters so that we are insensitive to where the data actually is
  VNetworkBase *GetNetAll() const;
  VNetworkBase *GetNetLow() const;
  VNetworkBase *GetNetHigh() const;

  bool LoadGANNoRange(std::string inputFileName);
  bool LoadGANFromRange(std::string inputFileName, std::string energyRange);

  TFCSGANXMLParameters m_param;

  ClassDef(TFCSGANEtaSlice, 5) // TFCSGANEtaSlice
};

#endif //> !ISF_TFCSGANETASLICE_H
