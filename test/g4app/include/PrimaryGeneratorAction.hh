// Copyright (c) 2024 CERN for the benefit of the FastCaloSim project

#ifndef PRIMARYGENERATORACTION_HH
#define PRIMARYGENERATORACTION_HH

#include <G4Types.hh>
#include "G4VUserPrimaryGeneratorAction.hh"
#include "TestHelpers/ParticleContainer.h"

class G4ParticleGun;
class G4ParticleDefinition;

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
    PrimaryGeneratorAction();
    virtual ~PrimaryGeneratorAction() override;

    virtual void GeneratePrimaries(G4Event *) override;

    void configure_gun(const TestHelpers::Particle* particle);
    void set_particle_container(const TestHelpers::ParticleContainer* container){fParticleContainer = container;};

private:
    G4ParticleGun *fParticleGun;
    const TestHelpers::ParticleContainer *fParticleContainer;
};

#endif
