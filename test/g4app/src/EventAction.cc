// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#include "G4Event.hh"

#include "EventAction.hh"

EventAction::EventAction() = default;

EventAction::~EventAction()
= default;

void EventAction::BeginOfEventAction(const G4Event * /*anEvent*/)
{
}

void EventAction::EndOfEventAction(const G4Event *anEvent)
{
    auto *HCE = anEvent->GetHCofThisEvent();
    if (HCE == nullptr) {
        return;
}
}
