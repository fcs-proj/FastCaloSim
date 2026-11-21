// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#include "ActionInitialization.hh"

#include "EventAction.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"

ActionInitialization::ActionInitialization() = default;

ActionInitialization::~ActionInitialization() = default;

void ActionInitialization::BuildForMaster() const
{
  SetUserAction(new RunAction);
}

void ActionInitialization::Build() const
{
  auto* primaryGeneratorAction = new PrimaryGeneratorAction;

  if (fParticleContainer == nullptr) {
    G4Exception("ActionInitialization::Build()",
                "ActionInitialization",
                FatalException,
                "Particle container not set");
  }

  primaryGeneratorAction->set_particle_container(fParticleContainer);

  SetUserAction(primaryGeneratorAction);

  SetUserAction(new RunAction);
  SetUserAction(new EventAction);
}
