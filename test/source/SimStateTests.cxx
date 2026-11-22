// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#include <array>
#include <cstddef>
#include <cstdint>
#include <set>
#include <utility>
#include <vector>

#include "SimStateTests.h"

#include <CLHEP/Random/RanluxEngine.h>
#include <gtest/gtest.h>

#include "FastCaloSim/Core/TFCSSimulationState.h"
#include "FastCaloSim/Geometry/Cell.h"

TEST_F(TFCSSimulationStateTest, Initialization)
{
  TFCSSimulationState sim_state;
  EXPECT_EQ(sim_state.randomEngine(), nullptr);
  EXPECT_EQ(sim_state.E(), 0);
  EXPECT_FALSE(sim_state.Ebin().has_value());
}

TEST_F(TFCSSimulationStateTest, EnergyManipulations)
{
  TFCSSimulationState sim_state;
  sim_state.set_Ebin(3);
  EXPECT_EQ(sim_state.Ebin(), 3U);
  sim_state.set_E(1, 100.0);
  EXPECT_DOUBLE_EQ(sim_state.E(1), 100.0);
  sim_state.add_E(1, 50.0);
  EXPECT_DOUBLE_EQ(sim_state.E(1), 150.0);
  sim_state.set_Efrac(1, 0.75);
  EXPECT_DOUBLE_EQ(sim_state.Efrac(1), 0.75);
}

TEST_F(TFCSSimulationStateTest, LayerEnergies)
{
  TFCSSimulationState sim_state;
  // Layers that were not set read as 0, also beyond the highest set layer
  EXPECT_EQ(sim_state.E(0), 0);
  EXPECT_EQ(sim_state.Efrac(23), 0);
  sim_state.set_E(12, 100.0);
  sim_state.set_Efrac(5, 0.5);
  EXPECT_DOUBLE_EQ(sim_state.E(12), 100.0);
  EXPECT_DOUBLE_EQ(sim_state.Efrac(5), 0.5);
  for (std::size_t layer = 0; layer < 24; ++layer) {
    if (layer != 12) {
      EXPECT_EQ(sim_state.E(layer), 0) << "layer " << layer;
    }
    if (layer != 5) {
      EXPECT_EQ(sim_state.Efrac(layer), 0) << "layer " << layer;
    }
  }
  // add_E to a layer that was not set starts from 0 and adds to the total
  sim_state.set_E(10.0);
  sim_state.add_E(3, 25.0);
  EXPECT_DOUBLE_EQ(sim_state.E(3), 25.0);
  EXPECT_DOUBLE_EQ(sim_state.E(), 35.0);
  // A lower layer set later keeps the higher ones
  sim_state.set_E(1, 7.0);
  EXPECT_DOUBLE_EQ(sim_state.E(1), 7.0);
  EXPECT_DOUBLE_EQ(sim_state.E(12), 100.0);
}

TEST_F(TFCSSimulationStateTest, ClearFunction)
{
  CLHEP::RanluxEngine engine(1);
  TFCSSimulationState sim_state(&engine);
  sim_state.set_Ebin(2);
  sim_state.set_E(100.0);
  sim_state.set_E(0, 100.0);
  sim_state.set_Efrac(0, 1.0);
  sim_state.deposit(7, 50.0);
  sim_state.set_SF(2.0);
  sim_state.setAuxInfo<int>("EventNr"_FCShash, 3);
  sim_state.setAuxObject<int>("BSEventData"_FCShash, 4);

  sim_state.clear();

  EXPECT_FALSE(sim_state.Ebin().has_value());
  EXPECT_EQ(sim_state.E(), 0);
  EXPECT_EQ(sim_state.E(0), 0);
  EXPECT_EQ(sim_state.Efrac(0), 0);
  EXPECT_TRUE(sim_state.cells().empty());
  EXPECT_EQ(sim_state.get_SF(), 1);
  EXPECT_FALSE(sim_state.hasAuxInfo("EventNr"_FCShash));
  EXPECT_EQ(sim_state.getAuxObject<int>("BSEventData"_FCShash), nullptr);
  EXPECT_EQ(sim_state.randomEngine(), &engine);
}

TEST_F(TFCSSimulationStateTest, DepositEnergy)
{
  TFCSSimulationState sim_state;
  // Define a new dummy cell
  Cell const cell(0,
                  Position {},
                  0,
                  /*isBarrel=*/true,
                  /*isXYZ=*/false,
                  /*isEtaPhiR=*/true,
                  /*isEtaPhiZ=*/false,
                  /*isRPhiZ=*/false,
                  0,
                  0,
                  0,
                  0,
                  0,
                  0);

  // Deposit 10 MeV in that cell
  sim_state.deposit(cell.id(), 10.0);
  // Retrieve the cell map of the simulation state
  const auto& cells = sim_state.cells();

  // Check that the cell is in the map
  ASSERT_TRUE(cells.contains(cell.id()));
  // Check that the energy deposited in the cell is correct
  EXPECT_DOUBLE_EQ(cells.at(cell.id()), 10.0);
}

TEST_F(TFCSSimulationStateTest, DepositSequences)
{
  TFCSSimulationState sim_state;
  // The same cell several times in a row, then another one and back
  sim_state.deposit(7, 1.0);
  sim_state.deposit(7, 2.0);
  sim_state.deposit(8, 4.0);
  sim_state.deposit(7, 8.0);
  EXPECT_DOUBLE_EQ(sim_state.cells().at(7), 11.0);
  EXPECT_DOUBLE_EQ(sim_state.cells().at(8), 4.0);

  // Many new cells (the map rehashes) between deposits into the same cell
  for (unsigned long long id = 100; id < 1100; ++id) {
    sim_state.deposit(8, 1.0);
    sim_state.deposit(id, 1.0);
  }
  sim_state.deposit(1099, 1.0);
  EXPECT_DOUBLE_EQ(sim_state.cells().at(8), 1004.0);
  EXPECT_DOUBLE_EQ(sim_state.cells().at(1099), 2.0);
  EXPECT_EQ(sim_state.cells().size(), 1002U);

  // A copy has its own cells
  TFCSSimulationState copy = sim_state;
  copy.deposit(1099, 1.0);
  sim_state.deposit(1099, 4.0);
  EXPECT_DOUBLE_EQ(copy.cells().at(1099), 3.0);
  EXPECT_DOUBLE_EQ(sim_state.cells().at(1099), 6.0);

  // So does a moved state
  TFCSSimulationState moved = std::move(copy);
  moved.deposit(1099, 1.0);
  EXPECT_DOUBLE_EQ(moved.cells().at(1099), 4.0);

  // Scaled cells keep collecting energy
  sim_state.scale_cells([](unsigned long long id) -> double
                        { return id == 1099 ? 0.5 : 1.0; });
  EXPECT_DOUBLE_EQ(sim_state.cells().at(1099), 3.0);
  EXPECT_DOUBLE_EQ(sim_state.cells().at(8), 1004.0);
  sim_state.deposit(1099, 1.0);
  EXPECT_DOUBLE_EQ(sim_state.cells().at(1099), 4.0);

  // So do cells after clear() (whether or not it erases them)
  sim_state.clear();
  const auto energy = [&sim_state](unsigned long long id) -> double
  {
    const auto it = sim_state.cells().find(id);
    return it == sim_state.cells().end() ? 0.0 : it->second;
  };
  const double before = energy(1099);
  sim_state.deposit(1099, 1.0);
  sim_state.deposit(1099, 1.0);
  EXPECT_DOUBLE_EQ(energy(1099), before + 2.0);
}

TEST_F(TFCSSimulationStateTest, DepositAfterMove)
{
  // The moved-from state must not write into the cells that the move took
  // along
  TFCSSimulationState source;
  source.deposit(7, 1.0);
  TFCSSimulationState constructed(std::move(source));
  // NOLINTNEXTLINE(bugprone-use-after-move,hicpp-invalid-access-moved)
  source.deposit(7, 2.0);
  EXPECT_DOUBLE_EQ(constructed.cells().at(7), 1.0);
  EXPECT_DOUBLE_EQ(source.cells().at(7), 2.0);

  // Move assignment
  TFCSSimulationState assigned;
  assigned.deposit(7, 8.0);
  assigned = std::move(constructed);
  // NOLINTNEXTLINE(bugprone-use-after-move,hicpp-invalid-access-moved)
  constructed.deposit(7, 4.0);
  assigned.deposit(7, 16.0);
  EXPECT_DOUBLE_EQ(assigned.cells().at(7), 17.0);
  EXPECT_DOUBLE_EQ(constructed.cells().at(7), 4.0);
}

TEST_F(TFCSSimulationStateTest, DepositAfterClear)
{
  // The cell energy after clear() is the energy deposited after it (clear()
  // erases the cells) or the sum (it keeps them), never a write through a
  // pointer to an erased cell
  TFCSSimulationState sim_state;
  sim_state.deposit(7, 1.0);
  sim_state.clear();
  const bool erased = !sim_state.cells().contains(7);
  sim_state.deposit(7, 2.0);
  ASSERT_TRUE(sim_state.cells().contains(7));
  EXPECT_DOUBLE_EQ(sim_state.cells().at(7), erased ? 2.0 : 3.0);
}

TEST_F(TFCSSimulationStateTest, AuxiliaryInfoHandling)
{
  TFCSSimulationState sim_state;
  uint32_t const hash = TFCSSimulationState::getAuxIndex("testInfo");
  sim_state.setAuxInfo<double>(hash, 123.456);
  EXPECT_DOUBLE_EQ(sim_state.getAuxInfo<double>(hash), 123.456);
}

TEST_F(TFCSSimulationStateTest, LayerAuxIndicesAreDistinct)
{
  // The per-layer entries of TFCSPredictExtrapWeights and
  // TFCSVoxelHistoLateralCovarianceFluctuations can be in the same simulation
  // state: their indices must differ from each other, from the bare layer
  // numbers and from the other named entries
  constexpr std::size_t n_layers = 24;
  const std::array layer_names {"ExtrapWeight"_FCShash,
                                "VoxelFluctWeights"_FCShash};
  const std::array other_names {"BSEinit"_FCShash,
                                "BSELayer"_FCShash,
                                "BSEventData"_FCShash,
                                "BSNHits"_FCShash,
                                "EventNr"_FCShash,
                                "FCSHitChainNHits"_FCShash,
                                "GANlayer"_FCShash,
                                "SF"_FCShash};

  std::set<std::uint32_t> indices(other_names.begin(), other_names.end());
  for (std::size_t layer = 0; layer < n_layers; ++layer) {
    indices.insert(static_cast<std::uint32_t>(layer));
  }
  for (const auto name : layer_names) {
    for (std::size_t layer = 0; layer < n_layers; ++layer) {
      indices.insert(TFCSSimulationState::layerAuxIndex(name, layer));
    }
  }
  EXPECT_EQ(indices.size(),
            other_names.size() + ((1 + layer_names.size()) * n_layers));

  // A value and an object for the same layer do not overwrite each other
  TFCSSimulationState sim_state;
  sim_state.setAuxInfo<double>(
      TFCSSimulationState::layerAuxIndex("ExtrapWeight"_FCShash, 2), 0.25);
  sim_state.setAuxObject(
      TFCSSimulationState::layerAuxIndex("VoxelFluctWeights"_FCShash, 2),
      std::vector<double> {1.0, 2.0});
  EXPECT_DOUBLE_EQ(
      sim_state.getAuxInfo<double>(
          TFCSSimulationState::layerAuxIndex("ExtrapWeight"_FCShash, 2)),
      0.25);
  const auto* weights = sim_state.getAuxObject<std::vector<double>>(
      TFCSSimulationState::layerAuxIndex("VoxelFluctWeights"_FCShash, 2));
  ASSERT_NE(weights, nullptr);
  EXPECT_EQ(weights->size(), 2U);
}

namespace
{
// Counts its live instances, to check who destroys the stored objects
struct Counted
{
  static inline int live = 0;
  int value = 0;
  explicit Counted(int v)
      : value(v)
  {
    ++live;
  }
  Counted(const Counted& other)
      : value(other.value)
  {
    ++live;
  }
  Counted(Counted&& other) noexcept
      : value(other.value)
  {
    ++live;
  }
  auto operator=(const Counted&) -> Counted& = default;
  auto operator=(Counted&&) noexcept -> Counted& = default;
  ~Counted() { --live; }
};
}  // namespace

TEST_F(TFCSSimulationStateTest, AuxObjectsAreOwned)
{
  Counted::live = 0;
  {
    TFCSSimulationState sim_state;
    const auto index = "TestObject"_FCShash;
    EXPECT_EQ(sim_state.getAuxObject<Counted>(index), nullptr);

    sim_state.setAuxObject(index, Counted(1));
    EXPECT_EQ(Counted::live, 1);
    const auto* first = sim_state.getAuxObject<Counted>(index);
    ASSERT_NE(first, nullptr);
    EXPECT_EQ(first->value, 1);
    // An object of another type is not returned
    EXPECT_EQ(sim_state.getAuxObject<int>(index), nullptr);

    // Replacing destroys the old object
    sim_state.setAuxObject(index, Counted(2));
    EXPECT_EQ(Counted::live, 1);
    const auto* second = sim_state.getAuxObject<Counted>(index);
    ASSERT_NE(second, nullptr);
    EXPECT_EQ(second->value, 2);

    // Erasing destroys it
    sim_state.eraseAuxObject(index);
    EXPECT_EQ(Counted::live, 0);
    EXPECT_EQ(sim_state.getAuxObject<Counted>(index), nullptr);

    // The state destroys what is still stored
    sim_state.setAuxObject(index, Counted(3));
    EXPECT_EQ(Counted::live, 1);
  }
  EXPECT_EQ(Counted::live, 0);
}
