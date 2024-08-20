#include "ActionInitialization.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "EventAction.hh"


ActionInitialization::ActionInitialization()
    : G4VUserActionInitialization()
{
}

ActionInitialization::~ActionInitialization()
{
}

void ActionInitialization::BuildForMaster() const
{
    SetUserAction(new RunAction);
}

void ActionInitialization::Build() const
{
    auto * primaryGeneratorAction = new PrimaryGeneratorAction;

    if (!fParticleContainer)
    {
        G4Exception("ActionInitialization::Build()", "ActionInitialization", FatalException, "Particle container not set");
    }

    primaryGeneratorAction->set_particle_container(fParticleContainer);

    SetUserAction(primaryGeneratorAction);

    SetUserAction(new RunAction);
    SetUserAction(new EventAction);
}
