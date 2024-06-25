#pragma once

#include <gtest/gtest.h>

#include "FastCaloSim/Core/MLogging.h"
#include "TestHelpers/Event.h"
#include "TestHelpers/ParticleContainer.h"
#include "TestHelpers/ParticleSampler.h"

class AtlasSimTestConfig : public ::testing::Environment
{
public:
  // Extrapolation debug level
  inline static const MSG::Level EXTRAPOLATION_MSG_LEVEL = MSG::Level::INFO;

  // The list of events (event = particle container) to be transported
  inline static const std::vector<TestHelpers::Event> EVENTS = []()
  {
    // Vector of events to be processed
    std::vector<TestHelpers::Event> events;

    // Set up the ID-Calo sampler
    TestHelpers::IDCaloBoundarySampler sampler;

    // Sample set of photons with specified eta
    TestHelpers::ParticleContainer photons =
        sampler.discreteEtaSample(22, 32768, {0.2});

    // Construct the event with the particle container and event label
    TestHelpers::Event evt(photons, std::string("Photons"));
    evt.set_print_string("Photons");

    // Add the event to the list of events
    events.emplace_back(evt);

    return events;
  }();

  // Method to set the output directory for the test
  inline static void set_output_dir(const std::string& dir)
  {
    output_dir = dir;
  }

  // Path to the serialized simulation data
  inline static auto simulation_output_path() -> std::string
  {
    return output_dir + "sim_state_data.json";
  }

private:
  inline static std::string output_dir;
};
