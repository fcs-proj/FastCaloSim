/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#ifndef ISF_FASTCALOSIMPARAMETRIZATION_CALOGEOMETRYLOOKUP_H
#define ISF_FASTCALOSIMPARAMETRIZATION_CALOGEOMETRYLOOKUP_H

#include <iostream>
#include <map>
#include <vector>

#include <TMath.h>

#include "FastCaloSim/Geometry/FSmap.h"
#include "FastCaloSim/Geometry/Identifier.h"
#include "FastCaloSim/Geometry/MeanAndRMS.h"

class CaloDetDescrElement;
class TCanvas;
class TGraphErrors;

typedef std::map<Identifier, const CaloDetDescrElement*> t_cellmap;
typedef std::map<double, const CaloDetDescrElement*> t_eta_cellmap;

class CaloGeometryLookup
{
public:
  CaloGeometryLookup(int ind = 0);
  virtual ~CaloGeometryLookup();

  bool IsCompatible(const CaloDetDescrElement* cell);
  void add(const CaloDetDescrElement* cell);
  t_cellmap::size_type size() const { return m_cells.size(); };
  int index() const { return m_index; };
  void set_index(int ind) { m_index = ind; };
  void post_process();
  bool has_overlap(CaloGeometryLookup* ref);
  void merge_into_ref(CaloGeometryLookup* ref);
  // void CalculateTransformation();

  float mineta() const { return m_mineta; };
  float maxeta() const { return m_maxeta; };
  float minphi() const { return m_minphi; };
  float maxphi() const { return m_maxphi; };

  float mineta_raw() const { return m_mineta_raw; };
  float maxeta_raw() const { return m_maxeta_raw; };
  float minphi_raw() const { return m_minphi_raw; };
  float maxphi_raw() const { return m_maxphi_raw; };

#ifdef USE_GPU
  // the geometry for GPU EDM are constructed with following functions
  float mineta_correction() const { return m_mineta_correction; };
  float maxeta_correction() const { return m_maxeta_correction; };
  float minphi_correction() const { return m_minphi_correction; };
  float maxphi_correction() const { return m_maxphi_correction; };
#endif

  float minx() const { return m_mineta; };
  float maxx() const { return m_maxeta; };
  float miny() const { return m_minphi; };
  float maxy() const { return m_maxphi; };

  float minx_raw() const { return m_mineta_raw; };
  float maxx_raw() const { return m_maxeta_raw; };
  float miny_raw() const { return m_minphi_raw; };
  float maxy_raw() const { return m_maxphi_raw; };

  const MeanAndRMS& deta() { return m_deta; };
  const MeanAndRMS& dphi() { return m_dphi; };
  float mindeta() { return m_mindeta; };
  float mindphi() { return m_mindphi; };
  const MeanAndRMS& dx() { return m_deta; };
  const MeanAndRMS& dy() { return m_dphi; };
  float mindx() { return m_mindeta; };
  float mindy() { return m_mindphi; };

  const MeanAndRMS& eta_correction() { return m_eta_correction; };
  const MeanAndRMS& phi_correction() { return m_phi_correction; };
  const MeanAndRMS& x_correction() { return m_eta_correction; };
  const MeanAndRMS& y_correction() { return m_phi_correction; };

  int cell_grid_eta() const { return m_cell_grid_eta; };
  int cell_grid_phi() const { return m_cell_grid_phi; };
  void set_xy_grid_adjustment_factor(float factor)
  {
    m_xy_grid_adjustment_factor = factor;
  };

  virtual const CaloDetDescrElement* getDDE(float eta,
                                            float phi,
                                            float* distance = 0,
                                            int* steps = 0);

#ifdef USE_GPU
  // used by FCS-GPU
  float xy_grid_adjustment_factor() const
  {
    return m_xy_grid_adjustment_factor;
  };
  float deta_double() const { return m_deta_double; };
  float dphi_double() const { return m_dphi_double; };
  const std::vector<std::vector<const CaloDetDescrElement*> >* cell_grid()
  {
    return &m_cell_grid;
  };
#endif

protected:
  float neta_double() { return (maxeta_raw() - mineta_raw()) / deta().mean(); };
  float nphi_double() { return (maxphi_raw() - minphi_raw()) / dphi().mean(); };
  Int_t neta() { return TMath::Nint(neta_double()); };
  Int_t nphi() { return TMath::Nint(nphi_double()); };

  // FCal is not sufficiently regular to use submaps with regular mapping
  float nx_double() { return (maxx_raw() - minx_raw()) / mindx(); };
  float ny_double() { return (maxy_raw() - miny_raw()) / mindy(); };
  Int_t nx() { return TMath::Nint(TMath::Ceil(nx_double())); };
  Int_t ny() { return TMath::Nint(TMath::Ceil(ny_double())); };

  float m_xy_grid_adjustment_factor;

  int raw_eta_position_to_index(float eta_raw) const
  {
    return TMath::Floor((eta_raw - mineta_raw()) / m_deta_double);
  };
  int raw_phi_position_to_index(float phi_raw) const
  {
    return TMath::Floor((phi_raw - minphi_raw()) / m_dphi_double);
  };
  bool index_range_adjust(int& ieta, int& iphi) const;
  float calculate_distance_eta_phi(const CaloDetDescrElement* DDE,
                                   float eta,
                                   float phi,
                                   float& dist_eta0,
                                   float& dist_phi0) const;

  int m_index;
  t_cellmap m_cells;
  std::vector<std::vector<const CaloDetDescrElement*> > m_cell_grid;
  int m_cell_grid_eta, m_cell_grid_phi;
  float m_mineta, m_maxeta, m_minphi, m_maxphi;
  float m_mineta_raw, m_maxeta_raw, m_minphi_raw, m_maxphi_raw;
  float m_mineta_correction, m_maxeta_correction, m_minphi_correction,
      m_maxphi_correction;

  MeanAndRMS m_deta, m_dphi, m_eta_correction, m_phi_correction;
  float m_mindeta, m_mindphi;
  float m_deta_double, m_dphi_double;
};
#endif
