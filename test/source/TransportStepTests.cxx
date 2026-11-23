// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#include <cmath>
#include <vector>

#include <gtest/gtest.h>

#include "FastCaloSim/Transport/G4CaloTransportTool.h"
#include "G4Box.hh"
#include "G4DynamicParticle.hh"
#include "G4FieldTrack.hh"
#include "G4Gamma.hh"
#include "G4GeometryManager.hh"
#include "G4LogicalVolume.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4Track.hh"

// A neutral particle is transported with steps of at most 10 m. Where the next
// boundary is further away, the navigator returns kInfinity: the particle
// must then move by the proposed step, not by kInfinity.
TEST(TransportStepTest, NeutralParticleCrossesLargeVolume)
{
  G4Material* vacuum =
      G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");

  // A world of 100 m with the envelope, where the transport stops, from
  // x = 25 m to 35 m
  auto* world_box = new G4Box("World", 50 * m, 50 * m, 50 * m);
  auto* world_log = new G4LogicalVolume(world_box, vacuum, "World");
  auto* world = new G4PVPlacement(
      nullptr, G4ThreeVector(), world_log, "World", nullptr, false, 0);
  auto* envelope_box = new G4Box("Envelope", 5 * m, 5 * m, 5 * m);
  auto* envelope_log = new G4LogicalVolume(envelope_box, vacuum, "Envelope");
  new G4PVPlacement(nullptr,
                    G4ThreeVector(30 * m, 0, 0),
                    envelope_log,
                    "Envelope",
                    world_log,
                    false,
                    0);
  G4GeometryManager::GetInstance()->CloseGeometry(true, false, world);

  G4CaloTransportTool tool(world);
  tool.setTransportLimitVolume("Envelope");

  // A photon from the origin towards the envelope
  auto* photon = new G4DynamicParticle(
      G4Gamma::Definition(), G4ThreeVector(1, 0, 0), 10 * GeV);
  const G4Track track(photon, 0, G4ThreeVector(0, 0, 0));

  const std::vector<G4FieldTrack> steps = tool.transport(track);

  // Start, 10 m, 20 m, envelope at 25 m
  ASSERT_EQ(steps.size(), 4U);
  for (const G4FieldTrack& step : steps) {
    EXPECT_TRUE(std::isfinite(step.GetPosition().x()));
  }
  EXPECT_NEAR(steps[1].GetPosition().x(), 10 * m, 1e-6);
  EXPECT_NEAR(steps[2].GetPosition().x(), 20 * m, 1e-6);
  EXPECT_NEAR(steps[3].GetPosition().x(), 25 * m, 1e-6);
  EXPECT_NEAR(steps[3].GetPosition().y(), 0, 1e-6);

  G4GeometryManager::GetInstance()->OpenGeometry(world);
}
