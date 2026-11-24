// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#pragma once

#include <string>
#include <tuple>
#include <vector>

#include <RtypesCore.h>
#include <hdf5.h>

#include "FastCaloSim/Core/TFCSBinnedShowerBase.h"
#include "FastCaloSim/Core/TFCSLateralShapeParametrizationHitBase.h"
#include "FastCaloSim/Core/TFCSMLCalorimeterSimulator.h"
#include "FastCaloSim/Core/TFCSParametrizationBinnedChain.h"
#include "FastCaloSim/Core/TFCSSimulationState.h"

class CaloGeo;

class TFCSBinnedShowerONNX : public TFCSBinnedShowerBase
{
public:
  explicit TFCSBinnedShowerONNX(const char* name = nullptr,
                                const char* title = nullptr);

  ~TFCSBinnedShowerONNX() override;
  TFCSBinnedShowerONNX(const TFCSBinnedShowerONNX&) = delete;
  TFCSBinnedShowerONNX(TFCSBinnedShowerONNX&&) = delete;
  auto operator=(const TFCSBinnedShowerONNX&) -> TFCSBinnedShowerONNX& = delete;
  auto operator=(TFCSBinnedShowerONNX&&) -> TFCSBinnedShowerONNX& = delete;

  auto is_match_Ekin_bin(std::size_t /*Ekin_bin*/) const -> bool override
  {
    return true;
  }

  // Loads the ONNX model for the calorimeter simulation.
  auto load_simulator(const std::string& filename) -> bool
  {
    if (m_ai_simulator != nullptr) {
      delete m_ai_simulator;
      m_ai_simulator = nullptr;
    }
    m_ai_simulator = new TFCSMLCalorimeterSimulator();
    if (!m_ai_simulator->loadSimulator(filename)) {
      delete m_ai_simulator;
      m_ai_simulator = nullptr;
      return false;
    }
    return true;
  }

  // Loads the voxel boundaries and (potentially) the average showers for the
  // upscaling
  // TODO: Define HDF5 file format somewhere
  void load_meta_data(const std::string& filename,
                      std::vector<std::size_t>& layers);

protected:
  // Returns the event to be used
  auto get_event(TFCSSimulationState& simulstate,
                 double eta_center,
                 double phi_center,
                 double e_init,
                 std::size_t reference_layer_index) const
      -> FCSReturnCode override;

private:
  // Store a reference to the class that handles the actual ONNX call.
  TFCSMLCalorimeterSimulator* m_ai_simulator = nullptr;

  ClassDefOverride(TFCSBinnedShowerONNX, 2)  // TFCSBinnedShowerONNX
};
