// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#include "FastCaloSim/Core/TFCSBinnedShowerONNX.h"

#include "FastCaloSim/Core/TFCSBinnedShowerBase.h"
#include "FastCaloSim/Core/TFCSMLCalorimeterSimulator.h"
#include "FastCaloSim/Core/TFCSSimulationState.h"
#include "TBuffer.h"
#include "TClass.h"

//=============================================
//======= TFCSBinnedShowerONNX =========
//=============================================

TFCSBinnedShowerONNX::TFCSBinnedShowerONNX(const char* name, const char* title)
    : TFCSBinnedShowerBase(name, title)
{
}

TFCSBinnedShowerONNX::~TFCSBinnedShowerONNX()
{
  // Delete the AI simulator
  if (m_ai_simulator != nullptr) {
    delete m_ai_simulator;
    m_ai_simulator = nullptr;
  }
}

auto TFCSBinnedShowerONNX::get_event(TFCSSimulationState& simulstate,
                                     double eta_center,
                                     double phi_center,
                                     double e_init,
                                     std::size_t reference_layer_index) const
    -> FCSReturnCode
{
  (void)reference_layer_index;  // Unused parameter
  (void)phi_center;  // Unused parameter

  if (m_ai_simulator == nullptr) {
    FCS_MSG_ERROR(
        "No simulator loaded. Please load a model with load_simulator for "
        "TFCSBinnedShowerONNX::get_event.");
    return FCSFatal;
  }

  TFCSMLCalorimeterSimulator::event_t event =
      // The network takes its inputs in single precision
      m_ai_simulator->getEvent(simulstate,
                               static_cast<float>(eta_center),
                               static_cast<float>(e_init));

  store_event(simulstate, std::move(event.event_data), e_init);
  return FCSSuccess;
}

void TFCSBinnedShowerONNX::load_meta_data(const std::string& filename,
                                          std::vector<std::size_t>& layers)
{
  m_coordinates.clear();

  if (m_use_upscaling) {
    m_sub_bin_distribution.clear();
  }

  // layer dependent variables
  for (std::size_t const layer_index : layers) {
    FCS_MSG_INFO("Loading layer " << layer_index << " from file: " << filename);

    // Load the bin boundaries for this layer
    load_bin_boundaries(filename, layer_index);
  }
}

// TODO: Could probably be removed in the end
void TFCSBinnedShowerONNX::Streamer(TBuffer& buffer)
{
  // Stream an object of class TFCSBinnedShowerONNX

  if (buffer.IsReading()) {
    buffer.ReadClassBuffer(TFCSBinnedShowerONNX::Class(), this);

  } else {
    buffer.WriteClassBuffer(TFCSBinnedShowerONNX::Class(), this);
  }
}
