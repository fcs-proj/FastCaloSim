// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#include "RunAction.hh"

#include "G4RunManager.hh"

RunAction::RunAction() = default;

RunAction::~RunAction() = default;

void RunAction::BeginOfRunAction(const G4Run* /*unused*/) {}

void RunAction::EndOfRunAction(const G4Run* /*unused*/) {}
