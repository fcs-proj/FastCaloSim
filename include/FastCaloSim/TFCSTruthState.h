/*
  Copyright (C) 2002-2023 CERN for the benefit of the ATLAS collaboration
*/

#ifndef ISF_FASTCALOSIMEVENT_TFCSTruthState_h
#define ISF_FASTCALOSIMEVENT_TFCSTruthState_h


#include <FastCaloSim/FastCaloSim_export.h>

#include <TLorentzVector.h>

#include "FastCaloSim/MLogging.h"

class FASTCALOSIM_EXPORT TFCSTruthState
    : public TLorentzVector
    , public ISF_FCS::MLogging
{
public:
  TFCSTruthState();
  TFCSTruthState(Double_t x, Double_t y, Double_t z, Double_t t, int pdgid);

  void set_pdgid(int val) { m_pdgid = val; };
  void set_vertex(const TLorentzVector& val) { m_vertex = val; };
  void set_vertex(Double_t x, Double_t y, Double_t z, Double_t t = 0)
  {
    m_vertex.SetXYZT(x, y, z, t);
  };
  void set_Ekin_off(double val) { m_ekin_off = val; };

  int pdgid() const { return m_pdgid; };
  double Ekin() const { return E() - M() + m_ekin_off; };
  double Ekin_off() const { return m_ekin_off; };
  const TLorentzVector& vertex() const { return m_vertex; };

  void Print(Option_t* option = "") const;

private:
  int m_pdgid;
  TLorentzVector m_vertex;
  double m_ekin_off = 0;

  ClassDef(TFCSTruthState, 3)  // TFCSTruthState
};

#endif
