// Copyright (c) 2024 CERN for the benefit of the FastCaloSim project

#include <iostream>

#include "FastCaloSim/Core/TFCSTruthState.h"

//=============================================
//======= TFCSTruthState =========
//=============================================
using FastCaloSim::Core::TFCSTruthState;

TFCSTruthState::TFCSTruthState()
    : TLorentzVector()
    , m_pdgid(0)
    , m_vertex(0, 0, 0, 0)
{
}

TFCSTruthState::TFCSTruthState(
    Double_t x, Double_t y, Double_t z, Double_t t, int pdgid)
    : TLorentzVector(x, y, z, t)
    , m_vertex(0, 0, 0, 0)
{
  m_pdgid = pdgid;
}

void TFCSTruthState::Print(Option_t*) const
{
  MSG_INFO("PDGID=" << m_pdgid << " pT=" << Pt() << " eta=" << Eta() << " phi="
                    << Phi() << " E=" << E() << " Ekin_off=" << Ekin_off());
}
