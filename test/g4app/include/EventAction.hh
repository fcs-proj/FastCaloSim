#ifndef EVENTACTION_HH
#define EVENTACTION_HH

#include "G4UserEventAction.hh"

class EventAction : public G4UserEventAction
{
public:
    EventAction();
    ~EventAction() override;

    virtual void BeginOfEventAction(const G4Event *) override;
    virtual void EndOfEventAction(const G4Event *) override;

private:
};

#endif // EVENTACTION_HH
