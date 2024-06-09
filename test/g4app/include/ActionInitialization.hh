#ifndef ACTIONINITIALIZATION_HH
#define ACTIONINITIALIZATION_HH

#include <G4PrimaryParticle.hh>
#include "G4VUserActionInitialization.hh"
#include "TestHelpers/ParticleContainer.h"

class G4ParticleDefinition;

class ActionInitialization : public G4VUserActionInitialization
{

    
public:
    ActionInitialization();
    virtual ~ActionInitialization() override;

    virtual void BuildForMaster() const override;
    virtual void Build() const override;

    // @brief Set the container of particles you want to use
    void set_particle_container(const TestHelpers::ParticleContainer* container){fParticleContainer = container;};

    private:
    
    const TestHelpers::ParticleContainer *fParticleContainer;

};

#endif
