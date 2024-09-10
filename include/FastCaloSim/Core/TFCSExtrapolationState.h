/*
  Copyright (C) 2002-2023 CERN for the benefit of the ATLAS collaboration
*/

#ifndef ISF_FASTCALOSIMEVENT_TFCSExtrapolationState_h
#define ISF_FASTCALOSIMEVENT_TFCSExtrapolationState_h

#include <FastCaloSim/FastCaloSim_export.h>
#include <TObject.h>

#include "FastCaloSim/Core/MLogging.h"
#include "FastCaloSim/Geometry/FastCaloSim_CaloCell_ID.h"

class FASTCALOSIM_EXPORT TFCSExtrapolationState
    : public TObject
    , public ISF_FCS::MLogging
{
public:
  TFCSExtrapolationState();

  void clear();

  enum SUBPOS
  {
    SUBPOS_MID = 0,
    SUBPOS_ENT = 1,
    SUBPOS_EXT = 2
  };  // MID=middle, ENT=entrance, EXT=exit of cal layer

  void set_OK(int layer, int subpos, bool val = true)
  {
    m_CaloOK[layer][subpos] = val;
  };

  void set_eta(int layer, int subpos, double val)
  {
    m_etaCalo[layer][subpos] = val;
  };
  void set_phi(int layer, int subpos, double val)
  {
    m_phiCalo[layer][subpos] = val;
  };
  void set_r(int layer, int subpos, double val)
  {
    m_rCalo[layer][subpos] = val;
  };
  void set_z(int layer, int subpos, double val)
  {
    m_zCalo[layer][subpos] = val;
  };
  void set_d(int layer, int subpos, double val)
  {
    m_dCalo[layer][subpos] = val;
  };
  void set_detaBorder(int layer, int subpos, double val)
  {
    m_distetaCaloBorder[layer][subpos] = val;
  };

  void set_IDCaloBoundary_eta(double val) { m_IDCaloBoundary_eta = val; };
  void set_IDCaloBoundary_phi(double val) { m_IDCaloBoundary_phi = val; };
  void set_IDCaloBoundary_r(double val) { m_IDCaloBoundary_r = val; };
  void set_IDCaloBoundary_x(double val) { m_IDCaloBoundary_x = val; };
  void set_IDCaloBoundary_y(double val) { m_IDCaloBoundary_y = val; };
  void set_IDCaloBoundary_z(double val) { m_IDCaloBoundary_z = val; };

  bool OK(int layer, int subpos) const { return m_CaloOK[layer][subpos]; };
  double eta(int layer, int subpos) const { return m_etaCalo[layer][subpos]; };
  double phi(int layer, int subpos) const { return m_phiCalo[layer][subpos]; };
  double r(int layer, int subpos) const { return m_rCalo[layer][subpos]; };
  double z(int layer, int subpos) const { return m_zCalo[layer][subpos]; };
  double d(int layer, int subpos) const { return m_dCalo[layer][subpos]; };
  double detaBorder(int layer, int subpos) const
  {
    return m_distetaCaloBorder[layer][subpos];
  };

  double IDCaloBoundary_eta() const { return m_IDCaloBoundary_eta; };
  double IDCaloBoundary_phi() const { return m_IDCaloBoundary_phi; };
  double IDCaloBoundary_r() const { return m_IDCaloBoundary_r; };
  double IDCaloBoundary_x() const { return m_IDCaloBoundary_x; };
  double IDCaloBoundary_y() const { return m_IDCaloBoundary_y; };
  double IDCaloBoundary_z() const { return m_IDCaloBoundary_z; };

  double IDCaloBoundary_AngleEta() const { return m_IDCaloBoundary_AngleEta; };
  double IDCaloBoundary_Angle3D() const { return m_IDCaloBoundary_Angle3D; };

  void set_IDCaloBoundary_AngleEta(double val)
  {
    m_IDCaloBoundary_AngleEta = val;
  };
  void set_IDCaloBoundary_Angle3D(double val)
  {
    m_IDCaloBoundary_Angle3D = val;
  };

  void Print(Option_t* option = "") const;

private:
  bool m_CaloOK[CaloCell_ID_FCS::MaxSample][3];
  double m_etaCalo[CaloCell_ID_FCS::MaxSample][3];
  double m_phiCalo[CaloCell_ID_FCS::MaxSample][3];
  double m_rCalo[CaloCell_ID_FCS::MaxSample][3];
  double m_zCalo[CaloCell_ID_FCS::MaxSample][3];
  double m_dCalo[CaloCell_ID_FCS::MaxSample][3];
  double m_distetaCaloBorder[CaloCell_ID_FCS::MaxSample][3];

  double m_IDCaloBoundary_eta;
  double m_IDCaloBoundary_phi;
  double m_IDCaloBoundary_r;
  double m_IDCaloBoundary_x;
  double m_IDCaloBoundary_y;
  double m_IDCaloBoundary_z;

  double m_IDCaloBoundary_AngleEta;
  double m_IDCaloBoundary_Angle3D;

  ClassDef(TFCSExtrapolationState, 2)  // TFCSExtrapolationState
};

#endif
