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

    if (!fPrimaryParticle)
    {
        G4ExceptionDescription description;
        description << "Primary particle not set";
        G4Exception("ActionInitialization::Build()", "ActionInitialization", FatalException, description);
    }

    primaryGeneratorAction->set_particle(fPrimaryParticle);

    SetUserAction(primaryGeneratorAction);

    SetUserAction(new RunAction);
    SetUserAction(new EventAction);
}