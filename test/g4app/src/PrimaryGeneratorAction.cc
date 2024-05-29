#include "G4ParticleGun.hh"

#include "PrimaryGeneratorAction.hh"
#include "G4ParticleTable.hh"


PrimaryGeneratorAction::PrimaryGeneratorAction()
    : G4VUserPrimaryGeneratorAction()
{
  G4int n_particle = 1;
  fParticleGun = new G4ParticleGun(n_particle);
}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
    delete fParticleGun;
}

void PrimaryGeneratorAction::set_particle(const TestHelpers::Particle* particle)
{
    fParticleGun->SetParticleDefinition(G4ParticleTable::GetParticleTable()->FindParticle(particle->pid));
    fParticleGun->SetParticleEnergy(particle->ekin);
    fParticleGun->SetParticlePosition(particle->pos);
    fParticleGun->SetParticleMomentumDirection(particle->dir);
}
void PrimaryGeneratorAction::GeneratePrimaries(G4Event *anEvent)
{
    fParticleGun->GeneratePrimaryVertex(anEvent);
}