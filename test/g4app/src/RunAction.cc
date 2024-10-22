// Copyright (c) 2024 CERN for the benefit of the FastCaloSim project

#include "G4RunManager.hh"

#include "RunAction.hh"

RunAction::RunAction()
    : G4UserRunAction()
{
}

RunAction::~RunAction()
{
}

void RunAction::BeginOfRunAction(const G4Run *)
{
}

void RunAction::EndOfRunAction(const G4Run *)
{
}
