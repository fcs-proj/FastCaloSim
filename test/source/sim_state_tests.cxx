#include <CLHEP/Random/RanluxEngine.h>
#include <gtest/gtest.h>

#include "FastCaloSim/Core/TFCSSimulationState.h"
#include "FastCaloSim/Geometry/CaloDetDescrElement.h"

class TFCSSimulationStateTest : public ::testing::Test
{
protected:
  CLHEP::HepRandomEngine* m_random_engine;

  void SetUp() override { m_random_engine = new CLHEP::RanluxEngine(); }

  void TearDown() override { delete m_random_engine; }
};

TEST_F(TFCSSimulationStateTest, Initialization)
{
  TFCSSimulationState sim_state;
  EXPECT_EQ(sim_state.randomEngine(), nullptr);
  EXPECT_EQ(sim_state.E(), 0);
  EXPECT_EQ(sim_state.Ebin(), -1);
}

TEST_F(TFCSSimulationStateTest, EnergyManipulations)
{
  TFCSSimulationState sim_state;
  sim_state.set_Ebin(3);
  EXPECT_EQ(sim_state.Ebin(), 3);
  sim_state.set_E(1, 100.0);
  EXPECT_DOUBLE_EQ(sim_state.E(1), 100.0);
  sim_state.add_E(1, 50.0);
  EXPECT_DOUBLE_EQ(sim_state.E(1), 150.0);
  sim_state.set_Efrac(1, 0.75);
  EXPECT_DOUBLE_EQ(sim_state.Efrac(1), 0.75);
}

TEST_F(TFCSSimulationStateTest, ClearFunction)
{
  TFCSSimulationState sim_state;
  sim_state.set_Ebin(2);
  sim_state.set_E(0, 100.0);
  sim_state.clear();

  EXPECT_EQ(sim_state.Ebin(), -1);
  EXPECT_EQ(sim_state.E(), 0);
  for (int i = 0; i < CaloCell_ID_FCS::MaxSample; ++i) {
    EXPECT_EQ(sim_state.E(i), 0);
    EXPECT_EQ(sim_state.Efrac(i), 0);
  }
}

TEST_F(TFCSSimulationStateTest, DepositEnergy)
{
  TFCSSimulationState sim_state;
  // Define a new cell
  std::unique_ptr<CaloDetDescrElement> element(new CaloDetDescrElement());
  // Deposit 10 MeV in that cell
  sim_state.deposit(element.get(), 10.0);
  // Retrieve the cell map of the simulation state
  auto& cells = sim_state.cells();
  // Check that the cell is in the map
  ASSERT_TRUE(cells.find(element.get()) != cells.end());
  // Check that the energy deposited in the cell is correct
  EXPECT_FLOAT_EQ(cells[element.get()], 10.0);
}

TEST_F(TFCSSimulationStateTest, AuxiliaryInfoHandling)
{
  TFCSSimulationState sim_state;
  uint32_t hash = TFCSSimulationState::getAuxIndex("testInfo");
  sim_state.setAuxInfo<double>(hash, 123.456);
  EXPECT_DOUBLE_EQ(sim_state.getAuxInfo<double>(hash), 123.456);
}

auto main(int argc, char** argv) -> int
{
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
