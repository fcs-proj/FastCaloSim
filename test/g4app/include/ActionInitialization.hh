#ifndef ACTIONINITIALIZATION_HH
#define ACTIONINITIALIZATION_HH

#include <G4PrimaryParticle.hh>
#include "G4VUserActionInitialization.hh"
#include "TestHelpers/Particle.h"

class G4ParticleDefinition;

class ActionInitialization : public G4VUserActionInitialization
{

    
public:
    ActionInitialization();
    virtual ~ActionInitialization() override;

    virtual void BuildForMaster() const override;
    virtual void Build() const override;
    
    void set_particle(const TestHelpers::Particle* particle){fPrimaryParticle = particle;};
    const TestHelpers::Particle *fPrimaryParticle;

    private:
    
};

#endif
