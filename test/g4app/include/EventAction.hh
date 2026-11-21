// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#ifndef EVENTACTION_HH
#define EVENTACTION_HH

#include "G4UserEventAction.hh"

class EventAction : public G4UserEventAction
{
public:
  EventAction();
  ~EventAction() override;
  EventAction(const EventAction&) = delete;
  auto operator=(const EventAction&) -> EventAction& = delete;
  EventAction(EventAction&&) = delete;
  auto operator=(EventAction&&) -> EventAction& = delete;

  void BeginOfEventAction(const G4Event* /*anEvent*/) override;
  void EndOfEventAction(const G4Event* /*anEvent*/) override;

private:
};

#endif  // EVENTACTION_HH
