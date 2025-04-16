// Copyright (c) 2024 CERN for the benefit of the FastCaloSim project

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
}

TEST_F(TFCSSimulationStateTest, DepositEnergy)
{
  TFCSSimulationState sim_state;
  // Define a new dummy cell
  Cell cell(0,
            Position {0, 0, 0, 0, 0, 0},
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
  auto& cells = sim_state.cells();

  // Check that the cell is in the map
  ASSERT_TRUE(cells.find(cell.id()) != cells.end());
  // Check that the energy deposited in the cell is correct
  EXPECT_FLOAT_EQ(cells[cell.id()], 10.0);
}

TEST_F(TFCSSimulationStateTest, AuxiliaryInfoHandling)
{
  TFCSSimulationState sim_state;
  uint32_t hash = TFCSSimulationState::getAuxIndex("testInfo");
  sim_state.setAuxInfo<double>(hash, 123.456);
  EXPECT_DOUBLE_EQ(sim_state.getAuxInfo<double>(hash), 123.456);
}
