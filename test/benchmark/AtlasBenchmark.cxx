// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

// Throughput benchmark of the full FastCaloSim chain (transport, extrapolation
// and simulation) on the ATLAS test setup. For each AtlasSimulation test
// configuration it shoots BENCH_N particles (default 50), spread over phi and
// slightly in eta, and reports the wall-clock time per configuration.
//
// Every transport step, extrapolation value and simulated cell energy is
// written as a hexadecimal float to BENCH_DUMP (if set), so that two builds can
// be checked for bit-identical output with a plain diff. The dump depends on
// the random-number sequence, so only compare full runs with the same BENCH_N.

#include <cerrno>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <map>
#include <memory>
#include <print>
#include <stdexcept>
#include <string>

#include <gtest/gtest.h>

#include "AtlasSimTests.h"
#include "FastSimModel.hh"
#include "G4GlobalFastSimulationManager.hh"
#include "TestConfig/AtlasSimTestsConfig.h"
#include "TestHelpers/Event.h"
#include "TestHelpers/ParticleSampler.h"

namespace
{
auto env_int(const char* name, int fallback) -> int
{
  const char* value = std::getenv(name);
  return value != nullptr ? static_cast<int>(std::strtol(value, nullptr, 10))
                          : fallback;
}

// Single dump file shared by all configurations, opened on first use; null if
// BENCH_DUMP is not set
auto dump_file() -> FILE*
{
  static std::unique_ptr<FILE, int (*)(FILE*)> const file =
      []() -> std::unique_ptr<FILE, int (*)(FILE*)>
  {
    const char* path = std::getenv("BENCH_DUMP");
    FILE* f = path != nullptr ? std::fopen(path, "w") : nullptr;
    if (path != nullptr && f == nullptr) {
      throw std::runtime_error(std::string("Cannot open BENCH_DUMP file ")
                               + path + ": " + std::strerror(errno));
    }
    return {f, &std::fclose};
  }();
  return file.get();
}

// The floats are written with printf's %a: its hexadecimal format (with the 0x
// prefix) is the one of the reference dumps, while std::print's {:a} omits the
// prefix. A failed write shows up when comparing the dumps.
void dump(FILE* out, const std::string& label, FastSimModel& model)
{
  std::println(out, "CONFIG {}", label);

  // The getters return copies: keep them alive while iterating
  const auto tracks = model.GetTransportTracks();
  const auto extrapolations = model.GetExtrapolations();
  const auto states = model.GetSimulationStates();

  int itrack = 0;
  for (const auto& track : tracks.get()) {
    std::println(out, "TRACK {}", itrack++);
    for (const auto& pos : track) {
      // NOLINTNEXTLINE(modernize-use-std-print,cert-err33-c): %a, see above
      std::fprintf(out, "%a %a %a\n", pos.x(), pos.y(), pos.z());
    }
  }

  int iextrap = 0;
  for (const auto& ex : extrapolations.get()) {
    // NOLINTNEXTLINE(modernize-use-std-print,cert-err33-c): %a, see above
    std::fprintf(out,
                 "EXTRAP %d %a %a %a %a %a %a %a %a\n",
                 iextrap++,
                 ex.IDCaloBoundary_eta(),
                 ex.IDCaloBoundary_phi(),
                 ex.IDCaloBoundary_r(),
                 ex.IDCaloBoundary_x(),
                 ex.IDCaloBoundary_y(),
                 ex.IDCaloBoundary_z(),
                 ex.IDCaloBoundary_AngleEta(),
                 ex.IDCaloBoundary_Angle3D());
    for (std::size_t layer = 0; layer <= 23; ++layer) {
      for (const FastCaloSim::SubPos subpos : {FastCaloSim::SubPos::MID,
                                               FastCaloSim::SubPos::ENT,
                                               FastCaloSim::SubPos::EXT})
      {
        // NOLINTNEXTLINE(modernize-use-std-print,cert-err33-c): %a, see above
        std::fprintf(out,
                     "%zu %d %d %a %a %a %a %a %a\n",
                     layer,
                     static_cast<int>(subpos),
                     static_cast<int>(ex.OK(layer, subpos)),
                     ex.eta(layer, subpos),
                     ex.phi(layer, subpos),
                     ex.r(layer, subpos),
                     ex.z(layer, subpos),
                     ex.d(layer, subpos),
                     ex.detaBorder(layer, subpos));
      }
    }
  }

  int istate = 0;
  for (const auto& state : states.get()) {
    // -1 for no energy bin, as in the reference dumps
    if (const auto ebin = state.Ebin()) {
      // NOLINTNEXTLINE(modernize-use-std-print,cert-err33-c): %a, see above
      std::fprintf(out, "STATE %d %zu %a\n", istate++, *ebin, state.E());
    } else {
      // NOLINTNEXTLINE(modernize-use-std-print,cert-err33-c): %a, see above
      std::fprintf(out, "STATE %d -1 %a\n", istate++, state.E());
    }
    for (std::size_t sample = 0; sample <= 23; ++sample) {
      // NOLINTNEXTLINE(modernize-use-std-print,cert-err33-c): %a, see above
      std::fprintf(
          out, "E %zu %a %a\n", sample, state.E(sample), state.Efrac(sample));
    }
    const std::map<unsigned long long, float> cells(state.cells().begin(),
                                                    state.cells().end());
    for (const auto& [id, energy] : cells) {
      // NOLINTNEXTLINE(modernize-use-std-print,cert-err33-c): %a, see above
      std::fprintf(out, "%llu %a\n", id, static_cast<double>(energy));
    }
  }
}
class AtlasBenchmark : public AtlasSimTests
{
};
}  // namespace

TEST_P(AtlasBenchmark, Throughput)
{
  const TestHelpers::Event& evt = GetParam();
  const int n_particles = env_int("BENCH_N", 50);
  ASSERT_GT(n_particles, 0) << "BENCH_N must be a positive integer";
  // Open the dump before the run, so that a bad path fails without waiting
  FILE* const out = dump_file();

  // Replicate the configuration's particle, rotated in phi and spread by
  // +-0.05 in eta, so that transport and cell lookups see varied positions
  const TestHelpers::Particle& ref = evt.get_container().get().front();
  TestHelpers::IDCaloBoundarySampler sampler;
  TestHelpers::ParticleContainer particles;
  for (int i = 0; i < n_particles; ++i) {
    const double frac = (i + 0.5) / n_particles;
    const double eta = ref.dir.eta() + 0.1 * (frac - 0.5);
    const double phi = 2 * M_PI * frac;
    TestHelpers::Particle particle =
        sampler.generate(ref.pid, ref.ekin, static_cast<float>(eta));
    particle.pos.rotateZ(phi);
    particle.dir.rotateZ(phi);
    particles.add(particle);
  }

  G4RunTests::actionInitialization->set_particle_container(&particles);
  G4RunTests::run_manager->SetUserInitialization(
      G4RunTests::actionInitialization);
  G4RunTests::run_manager->Initialize();

  auto* mgr = G4GlobalFastSimulationManager::GetGlobalFastSimulationManager();
  auto* model =
      dynamic_cast<FastSimModel*>(mgr->GetFastSimulationModel("model"));
  ASSERT_NE(model, nullptr);
  model->Reset();

  BasicExtrapolTests::extrapolator.setLevel(FastCaloSim::MSG::Level::INFO);
  model->setExtrapolationTool(BasicExtrapolTests::extrapolator);
  model->doSimulation(/*doSimulation=*/true);

  auto* param = static_cast<TFCSParametrizationBase*>(
      BasicSimTests::param_files[evt.get_param_file_key()]->Get("SelPDGID"));
  ASSERT_NE(param, nullptr);
  param->set_geometry(AtlasGeoTests::geo);
  param->setLevel(FastCaloSim::MSG::Level::INFO, /*recursive=*/true);
  model->setParametrization(param);

  // Run initialization (physics tables on the first run) outside the timing
  G4RunTests::run_manager->BeamOn(0);

  const auto start = std::chrono::steady_clock::now();
  G4RunTests::run_manager->BeamOn(1);
  const auto stop = std::chrono::steady_clock::now();

  const double ms =
      std::chrono::duration<double, std::milli>(stop - start).count();
  std::println("BENCH {:<45} N={} ms={:.1f} ms/particle={:.3f}",
               evt.get_print_string(),
               n_particles,
               ms,
               ms / n_particles);

  if (out != nullptr) {
    dump(out, evt.get_print_string(), *model);
  }
}

INSTANTIATE_TEST_SUITE_P(AtlasBenchmark,
                         AtlasBenchmark,
                         ::testing::ValuesIn(AtlasSimTestConfig::EVENTS),
                         testing::PrintToStringParamName());
