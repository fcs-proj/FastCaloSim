// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#ifndef TFCSMLCALORIMETERSIMULATOR_H
#define TFCSMLCALORIMETERSIMULATOR_H

#include <string>
#include <vector>

#include "FastCaloSim/Core/MLogging.h"
#include "FastCaloSim/Core/TFCSSimulationState.h"

// generic network class
#include "FastCaloSim/Core/VNetworkBase.h"

class TFCSMLCalorimeterSimulator : public ISF_FCS::MLogging
{
public:
  TFCSMLCalorimeterSimulator();
  virtual ~TFCSMLCalorimeterSimulator();

  typedef struct
  {
    std::vector<unsigned int> bin_index_vector;
    std::vector<float> E_vector;
  } layer_t;

  typedef struct
  {
    std::vector<layer_t> event_data;
  } event_t;

  bool loadSimulator(std::string filename);

  void Print() const;

  VNetworkBase::NetworkOutputs predictVoxels(TFCSSimulationState& simulstate,
                                             float eta,
                                             float energy) const;
  event_t getEvent(TFCSSimulationState& simulstate,
                   float eta,
                   float energy) const;
  VNetworkBase::NetworkOutputs predictVoxels() const;

  void setInputShapes(std::vector<long unsigned int> layer_boundaries,
                      std::vector<long unsigned int> used_layers)
  {
    if (layer_boundaries.empty() || used_layers.empty()) {
      FCS_MSG_ERROR(
          "setInputShapes requires non-empty layer_boundaries "
          "and used_layers");
      return;
    }
    if (layer_boundaries.size() != used_layers.size() + 1) {
      FCS_MSG_ERROR("setInputShapes: layer_boundaries size ("
                    << layer_boundaries.size()
                    << ") must equal used_layers size (" << used_layers.size()
                    << ") + 1");
      return;
    }
    m_layer_boundaries = std::move(layer_boundaries);
    m_used_layers = std::move(used_layers);
    m_nVoxels = m_layer_boundaries.back();
    m_nLayers = m_used_layers.size();
  };

private:
  std::unique_ptr<VNetworkBase> m_onnx_model = nullptr;

  int m_nEvents = 1;  // Currently no batching supported from ONNX handler

  // Default shapes for the photon barrel-CFM
  // Should be set using setInputShapes for other models
  std::vector<long unsigned int> m_layer_boundaries = {
      0, 36, 200, 310, 346, 382};
  std::vector<long unsigned int> m_used_layers = {0, 1, 2, 3, 12};
  long unsigned int m_nVoxels = 382;
  long unsigned int m_nLayers = 5;

  ClassDef(TFCSMLCalorimeterSimulator, 1)  // TFCSMLCalorimeterSimulator
};

#endif  //> !ISF_TFCSMLCALORIMETERSIMULATOR_H
