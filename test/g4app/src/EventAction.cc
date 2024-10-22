// Copyright (c) 2024 CERN for the benefit of the FastCaloSim project

#include "G4Event.hh"

#include "EventAction.hh"

EventAction::EventAction()
    : G4UserEventAction()
{
}

EventAction::~EventAction()
{
}

void EventAction::BeginOfEventAction(const G4Event *)
{
}

void EventAction::EndOfEventAction(const G4Event *anEvent)
{
    auto HCE = anEvent->GetHCofThisEvent();
    if (!HCE)
        return;
}
