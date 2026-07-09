// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#pragma once

#include <string>

#include "FastCaloSim/Core/TFCSLateralShapeParametrizationHitBase.h"
#include "FastCaloSim/Core/TFCSSimulationState.h"

// forward declare lwtnn dependencies
namespace lwt
{
class LightweightNeuralNetwork;
}

class TFCSPredictExtrapWeights : public TFCSLateralShapeParametrizationHitBase
{
public:
  TFCSPredictExtrapWeights(const char* name = nullptr,
                           const char* title = nullptr,
                           CaloGeo* geo = nullptr);

  void set_geometry(CaloGeo* geo) override { m_geo = geo; };

  ~TFCSPredictExtrapWeights() override;

  auto operator==(const TFCSParametrizationBase& ref) const -> bool override;

  // Used to decorate simulstate with extrapolation weights
  auto simulate(TFCSSimulationState& simulstate,
                const TFCSTruthState* truth,
                const TFCSExtrapolationState* extrapol) const
      -> FCSReturnCode override;

  // Used to decorate Hit with extrapolated center positions
  auto simulate_hit(Hit& hit,
                    TFCSSimulationState& simulstate,
                    const TFCSTruthState* truth,
                    const TFCSExtrapolationState* extrapol)
      -> FCSReturnCode override;

  // Status bit for chain persistency
  enum FCSfreemem
  {
    kfreemem = BIT(17)  ///< Set this bit in the TObject bit if the memory for
                        ///< m_input should be freed after reading in athena
  };
  auto freemem() const -> bool { return TestBit(kfreemem); };
  void set_freemem() { SetBit(kfreemem); };

  // Initialize Neural Network
  auto initializeNetwork(int pid,
                         const std::string& etaBin,
                         const std::string& FastCaloNNInputFolderName) -> bool;

  // Get inputs needed to normalize data
  auto getNormInputs(const std::string& etaBin,
                     const std::string& FastCaloTXTInputFolderName) -> bool;

  // Prepare inputs to the Neural Network
  auto prepareInputs(TFCSSimulationState& simulstate, const float truthE) const
      -> std::map<std::string, double>;

  // Print()
  void Print(Option_t* option = "") const override;

  // Use extrapWeight=0.5 or r and z when constructing a hit?
  enum TFCSPredictExtrapWeightsStatusBits
  {
    kUseHardcodedWeight = BIT(15)
  };
  auto UseHardcodedWeight() const -> bool
  {
    return TestBit(kUseHardcodedWeight);
  };
  void set_UseHardcodedWeight() { SetBit(kUseHardcodedWeight); };
  void reset_UseHardcodedWeight() { ResetBit(kUseHardcodedWeight); };

protected:
  CaloGeo* m_geo;  //! do not persistify

private:
  // Persistify configuration in string m_input. A custom Streamer(...) builds
  // m_nn on the fly when reading from file.
  // Inside Athena, if freemem() is true, the content of m_input is deleted
  // after reading in order to free memory
  std::string* m_input = nullptr;
  std::vector<int>* m_relevantLayers = nullptr;
  lwt::LightweightNeuralNetwork* m_nn = nullptr;  //! Do not persistify
  std::vector<int>* m_normLayers =
      nullptr;  // vector of index layers (-1 corresponds to truth energy)
  std::vector<float>* m_normMeans =
      nullptr;  // vector of mean values for normalizing energy fraction per
                // layer, last index is for total energy
  std::vector<float>* m_normStdDevs =
      nullptr;  // vector of std dev values for normalizing energy fraction per
                // layer, last index is for total energy

  ClassDefOverride(TFCSPredictExtrapWeights, 1)  // TFCSPredictExtrapWeights
};
