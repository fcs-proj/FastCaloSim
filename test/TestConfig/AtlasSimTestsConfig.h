// Copyright (c) 2024 CERN for the benefit of the FastCaloSim project

#pragma once

#include <gtest/gtest.h>

#include "FastCaloSim/Core/MLogging.h"
#include "TestHelpers/Event.h"
#include "TestHelpers/IOManager.h"
#include "TestHelpers/ParticleContainer.h"
#include "TestHelpers/ParticleSampler.h"
#include "TestHelpers/ParticleSimParams.h"

class AtlasSimTestConfig : public ::testing::Environment
{
public:
  // Extrapolation debug level
  inline static const MSG::Level EXTRAPOLATION_MSG_LEVEL = MSG::Level::INFO;

  // The list of events (event = particle container) to be simulated
  inline static const std::vector<TestHelpers::Event> EVENTS = []()
  {
    // Vector of events to be processed
    std::vector<TestHelpers::Event> events;

    // Set up the ID-Calo sampler
    TestHelpers::IDCaloBoundarySampler sampler;

    // Define the properties of the particles you want to test
    // Alias for a vector of ParticleParams
    using ParticleParamsList = std::vector<ParticleSimParams>;

    // Define the properties of the particles you want to test
    ParticleParamsList params = {
        // Barrel (eta = 0.2)
        {22, 16384, 0.2, "barrel"},
        {22, 131072, 0.2, "barrel"},
        {211, 16384, 0.2, "barrel"},
        {211, 131072, 0.2, "barrel"},
        // Endcap (eta = 2.0)
        {22, 16384, 2.0, "endcap"},
        {22, 131072, 2.0, "endcap"},
        {211, 16384, 2.0, "endcap"},
        {211, 131072, 2.0, "endcap"},
        // FCal (eta = 3.5)
        {22, 16384, 3.5, "fcal"},
        {22, 131072, 3.5, "fcal"},
        {211, 16384, 3.5, "fcal"},
        {211, 131072, 3.5, "fcal"},
        // Barrel-Endcap Transition (eta = 1.45)
        {22, 16384, 1.45, "barrel_endcap_transition"},
        {22, 131072, 1.45, "barrel_endcap_transition"},
        {211, 16384, 1.45, "barrel_endcap_transition"},
        {211, 131072, 1.45, "barrel_endcap_transition"},
        // Endcap-FCal Transition (eta = 3.2)
        {22, 16384, 3.2, "endcap_fcal_transition"},
        {22, 131072, 3.2, "endcap_fcal_transition"},
        {211, 16384, 3.2, "endcap_fcal_transition"},
        {211, 131072, 3.2, "endcap_fcal_transition"}};

    // Loop over the parameters and create the events
    for (const auto& ptcl : params) {
      // Create photon with speficfic PID, EKIN, and ETA
      TestHelpers::Particle particle =
          sampler.generate(ptcl.pid, ptcl.ekin, ptcl.eta);

      // Create a particle container and add the particle
      TestHelpers::ParticleContainer particles;
      particles.add(particle);

      // Construct the event with the particle container and event label
      TestHelpers::Event evt(particles, ptcl.latex_label);
      // Set the name of the param file used to simulate the particle
      evt.set_param_file_key(ptcl.param_file_key);
      evt.set_print_string(ptcl.print_label);

      // Add the event to the list of events
      events.emplace_back(evt);
    }

    return events;
  }();

  // The location of the python script to plot the colored cells containing
  // energy
  inline static const std::string PYTHON_SCRIPT =
      std::string(TEST_BASE_DIR) + "/python/scripts/plot_simulation.py";

  // Set the name of the output files
  inline static const std::string TRANSPORT_DATA_FILE_NAME =
      "transport_data.json";
  inline static const std::string SIM_DATA_FILE_NAME = "sim_state_data.json";
  inline static const std::string SIM_PLOT_FILE_NAME =
      "simulated_cell_energy.png";
  inline static const std::string SIM_GIF_FILE_NAME =
      "simulated_cell_energy.gif";

  // Path to the serialized transportation data
  inline static auto transport_output_path() -> std::string
  {
    return output_dir + TRANSPORT_DATA_FILE_NAME;
  }

  // Path to the reference transportation data
  inline static auto transport_output_ref_path() -> std::string
  {
    return TestHelpers::IOManager::ref_dir() + TRANSPORT_DATA_FILE_NAME;
  }

  // Path to the serialized simulation data
  inline static auto sim_output_path() -> std::string
  {
    return output_dir + SIM_DATA_FILE_NAME;
  }

  // Path to the reference simulation data
  inline static auto sim_output_ref_path() -> std::string
  {
    return TestHelpers::IOManager::ref_dir() + SIM_DATA_FILE_NAME;
  }

  // Path to the output plot of the simulation plot
  inline static auto plot_output_sim_path() -> std::string
  {
    return output_dir + SIM_PLOT_FILE_NAME;
  }

  // Path to the output simulation animation
  inline static auto plot_output_sim_animation() -> std::string
  {
    return output_dir + SIM_GIF_FILE_NAME;
  }

  // String of the python command to plot simulation output
  inline static auto py_sim_exec(std::string label) -> std::string
  {
    std::string command = "python3 " + PYTHON_SCRIPT + " --sim_input "
        + sim_output_path() + " --output_png " + plot_output_sim_path()
        + " --output_gif " + plot_output_sim_animation() + " --transport_input "
        + transport_output_path() + " --label \"" + label + "\"";

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
