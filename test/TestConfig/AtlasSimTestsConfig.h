#pragma once

#include <gtest/gtest.h>

#include "FastCaloSim/Core/MLogging.h"
#include "TestHelpers/Event.h"
#include "TestHelpers/ParticleContainer.h"
#include "TestHelpers/ParticleSampler.h"

enum PARAM
{
  PID,
  EKIN,
  ETA,
  LATEX_LABEL,
  PRINT_LABEL
};

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
    std::vector<std::tuple<int, int, double, std::string, std::string>> params =
        {std::make_tuple(22,
                         16384,
                         0.2,
                         R"($ E_{\gamma}=16\,\text{GeV},\, \eta=0.20$)",
                         "pid22_E16_eta20"),
         std::make_tuple(22,
                         65536,
                         0.2,
                         R"($ E_{\gamma}=65\,\text{GeV},\, \eta=0.20$)",
                         "pid22_E65_eta20"),
         std::make_tuple(22,
                         131072,
                         0.2,
                         R"($ E_{\gamma}=131\,\text{GeV},\, \eta=0.20$)",
                         "pid22_E131_eta20")};

    // Loop over the parameters and create the events
    for (const auto& param : params) {
      // Creat photon with speficfic PID, EKIN, and ETA
      TestHelpers::Particle particle = sampler.generate(
          std::get<PID>(param), std::get<EKIN>(param), std::get<ETA>(param));

      // Create a particle container and add the particle
      TestHelpers::ParticleContainer particles;
      particles.add(particle);

      // Construct the event with the particle container and event label
      TestHelpers::Event evt(particles, std::get<LATEX_LABEL>(param));
      evt.set_print_string(std::get<PRINT_LABEL>(param));

      // Add the event to the list of events
      events.emplace_back(evt);
    }

    return events;
  }();

  // The location of the python script to plot the coloured cells containing
  // energy
  inline static const std::string PYTHON_SCRIPT =
      std::string(TEST_BASE_DIR) + "/python/plot_simulation.py";

  // Path to the serialized simulation data
  inline static auto sim_output_path() -> std::string
  {
    return output_dir + "sim_state_data.json";
  }

  // Path to the output plot of the simulation plot
  inline static auto plot_output_sim_path() -> std::string
  {
    return output_dir + "simulated_cell_energy.png";
  }

  // String of the python command to plot simulation output
  inline static auto py_sim_exec(std::string label) -> std::string
  {
    std::string command = "python3 " + PYTHON_SCRIPT + " --input "
        + sim_output_path() + " --output " + plot_output_sim_path()
        + " --label \"" + label + "\"";

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
