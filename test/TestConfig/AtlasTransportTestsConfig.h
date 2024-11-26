// Copyright (c) 2024 CERN for the benefit of the FastCaloSim project

#pragma once

#include <gtest/gtest.h>

#include "FastCaloSim/Core/MLogging.h"
#include "TestHelpers/Event.h"
#include "TestHelpers/IOManager.h"
#include "TestHelpers/ParticleContainer.h"
#include "TestHelpers/ParticleSampler.h"
#include "TestHelpers/ParticleTypes.h"

class AtlasTransportTestConfig : public ::testing::Environment
{
public:
  // Boolean flag to zoom into region with tracks in plotting
  inline static const bool ZOOM_IN = false;

  // Extrapolation debug level
  inline static const FCS_MSG::Level EXTRAPOLATION_FCS_MSG_LEVEL =
      FCS_MSG::Level::INFO;

  // The output directory for the test
  inline static const std::string OUTPUT_DIR = "output/";

  // The output directory for the test
  inline static const std::string TEST_GEO_DIR = "test/geo/";

  // The output directory for the test
  inline static const std::string PYTHON_SCRIPT_DIR = "test/python/scripts/";

  // The list of events (event = particle container) to be transported
  inline static const std::vector<TestHelpers::Event> EVENTS = []()
  {
    // Vector of events to be processed
    std::vector<TestHelpers::Event> events;

    // For each event we sample particles of the same type uniformly in eta
    const float min_eta = -5;  //-5
    const float max_eta = 5;  // 5
    const float step_size = 0.05;  // 0.5
    const float energy = 50 * CLHEP::MeV;
    std::vector<TestHelpers::ParticleType> ptype = {
        TestHelpers::ParticleTypes::Photon,
        TestHelpers::ParticleTypes::Electron,
        TestHelpers::ParticleTypes::Positron,
        TestHelpers::ParticleTypes::Proton};

    // Particle sampler on the ID-calorimeter boundary
    TestHelpers::IDCaloBoundarySampler sampler;

    for (const auto& type : ptype) {
      TestHelpers::ParticleContainer particles = sampler.uniformEtaSample(
          type.pid, energy, min_eta, max_eta, step_size);

      // Plot label rounded to 2 decimal places
      std::stringstream estring;
      estring << std::fixed << std::setprecision(2) << energy;
      std::string evt_label = "$ E=" + estring.str() + R"(\,\text{MeV}\,)"
          + std::string(type.label) + "$";
      // Construct the event with the particle container and event label
      TestHelpers::Event evt(particles, evt_label);
      // How should the test for this event be printed?
      std::string ptcl_name = std::string(type.name);
      evt.set_print_string(ptcl_name);

      events.emplace_back(evt);
    }

    return events;
  }();

  // The location of the python script to plot the transport
  inline static const std::string PYTHON_SCRIPT =
      std::string(TEST_BASE_DIR) + "/python/scripts/plot_transport.py";

  // The location of the python script to plot the extrapolation
  inline static const std::string PYTHON_EXTRAP_SCRIPT =
      std::string(TEST_BASE_DIR) + "/python/scripts/plot_extrapolation.py";

  // The location of calo layer data used for plotting the transport
  inline static const std::string CALO_LAYER_DATA =
      std::string(TEST_BASE_DIR) + "/python/scripts/calo_layers.csv";

  // Set the name of the output files
  inline static const std::string TRANSPORT_DATA_FILE_NAME =
      "transported_tracks.json";
  inline static const std::string EXTRAPOLATION_DATA_FILE_NAME =
      "extrapolation.json";
  inline static const std::string TRANSPORT_PLOT_NAME =
      "particle_transport.png";
  inline static const std::string EXTRAPOLATION_PLOT_NAME =
      "particle_extrapolation.png";

  // Path to the serialized transport data
  inline static auto transport_output_path() -> std::string
  {
    return output_dir + TRANSPORT_DATA_FILE_NAME;
  }

  // Path to the reference transport data
  inline static auto transport_output_ref_path() -> std::string
  {
    return TestHelpers::IOManager::ref_dir() + TRANSPORT_DATA_FILE_NAME;
  }

  // Path to the serialized extrapolation data
  inline static auto extrapolation_output_path() -> std::string
  {
    return output_dir + EXTRAPOLATION_DATA_FILE_NAME;
  }

  // Path to the reference extrapolation data
  inline static auto extrapolation_output_ref_path() -> std::string
  {
    return TestHelpers::IOManager::ref_dir() + EXTRAPOLATION_DATA_FILE_NAME;
  }

  // Path to the output plot of the transport
  inline static auto plot_output_transport_path() -> std::string
  {
    return output_dir + TRANSPORT_PLOT_NAME;
  }

  // Path to the output plot of the extrapolation
  inline static auto plot_output_extrapol_path() -> std::string
  {
    return output_dir + EXTRAPOLATION_PLOT_NAME;
  }

  // String of the python command to plot transportation
  inline static auto py_transport_exec(std::string label) -> std::string
  {
    std::string command = "python3 " + PYTHON_SCRIPT + " --input "
        + transport_output_path() + " --output " + plot_output_transport_path()
        + " --calo_layer_csv " + CALO_LAYER_DATA + " --label \"" + label + "\"";

    if (ZOOM_IN) {
      command += " --track_zoom";
    }

    return command;
  }

  // String of the python command to plot the extrapolation
  inline static auto py_extrapol_exec(std::string label) -> std::string
  {
    std::string command = "python3 " + PYTHON_EXTRAP_SCRIPT + " --input "
        + extrapolation_output_path() + " --output "
        + plot_output_extrapol_path() + " --label \"" + label + "\"";

    return command;
  }

  // Method to set the output directory for the test
  inline static void set_output_dir(const std::string& dir)
  {
    output_dir = dir;
  }

private:
  inline static std::string output_dir;
};
