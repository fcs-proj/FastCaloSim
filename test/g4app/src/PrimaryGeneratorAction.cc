// Copyright (c) 2024 CERN for the benefit of the FastCaloSim project

#include "G4ParticleGun.hh"

#include "PrimaryGeneratorAction.hh"
#include "G4ParticleTable.hh"


PrimaryGeneratorAction::PrimaryGeneratorAction()
    : G4VUserPrimaryGeneratorAction()
{
  fParticleGun = new G4ParticleGun(1);
}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
    delete fParticleGun;
}

void PrimaryGeneratorAction::configure_gun(const TestHelpers::Particle* particle)
{
    fParticleGun->SetParticleDefinition(G4ParticleTable::GetParticleTable()->FindParticle(particle->pid));
    fParticleGun->SetParticleEnergy(particle->ekin);
    fParticleGun->SetParticlePosition(particle->pos);
    fParticleGun->SetParticleMomentumDirection(particle->dir);
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event *anEvent)
{

    if (fParticleContainer == nullptr)
    {
        G4Exception("PrimaryGeneratorAction::GeneratePrimaries()", "PrimaryGeneratorAction", FatalException, "Particle container not set");
    }

    for(const auto &particle : fParticleContainer->get())
    {
        configure_gun(&particle);
        fParticleGun->GeneratePrimaryVertex(anEvent);
    }
}
