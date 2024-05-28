#ifndef PRIMARYGENERATORACTION_HH
#define PRIMARYGENERATORACTION_HH

#include <G4Types.hh>
#include "G4VUserPrimaryGeneratorAction.hh"
#include "TestHelpers/Particle.h"

class G4ParticleGun;
class G4ParticleDefinition;

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
    PrimaryGeneratorAction();
    virtual ~PrimaryGeneratorAction() override;

    virtual void GeneratePrimaries(G4Event *) override;

    void set_particle(const TestHelpers::Particle* particle);

private:
    G4ParticleGun *fParticleGun;
};

#endif
