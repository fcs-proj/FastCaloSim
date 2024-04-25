/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#ifndef ISF_FASTCALOSIMPARAMETRIZATION_CALOGEOMETRY_H
#define ISF_FASTCALOSIMPARAMETRIZATION_CALOGEOMETRY_H

#include <atomic>
#include <iostream>
#include <map>
#include <vector>

#include <FastCaloSim/FastCaloSim_export.h>
#include <TMath.h>

#include "FastCaloSim/Geometry/CaloGeometryLookup.h"
#include "FastCaloSim/Geometry/FCAL_ChannelMap.h"
#include "FastCaloSim/Geometry/FSmap.h"
#include "FastCaloSim/Geometry/ICaloGeometry.h"
#include "FastCaloSim/Geometry/MeanAndRMS.h"

class CaloDetDescrElement;
class TCanvas;
class TGraph;
class TGraphErrors;

class FASTCALOSIM_EXPORT CaloGeometry : virtual public ICaloGeometry
{
public:
  static const int MAX_SAMPLING;

  static const Identifier m_debug_identify;
  static std::atomic<bool> m_debug;

  CaloGeometry();
  virtual ~CaloGeometry();

  virtual bool PostProcessGeometry();

  virtual void Validate(int nrnd = 100);

  virtual const CaloDetDescrElement* getDDE(Identifier identify);
  virtual const CaloDetDescrElement* getDDE(int sampling, Identifier identify);

  virtual const CaloDetDescrElement* getDDE(
      int sampling, float eta, float phi, float* distance = 0, int* steps = 0);
  virtual const CaloDetDescrElement* getFCalDDE(int sampling,
                                                float x,
                                                float y,
                                                float z,
                                                float* distance = 0,
                                                int* steps = 0);
  bool getClosestFCalCellIndex(
      int sampling, float x, float y, int& ieta, int& iphi, int* steps = 0);

  double deta(int sample, double eta) const;
  void minmaxeta(int sample, double eta, double& mineta, double& maxeta) const;
  double rzmid(int sample, double eta) const;
  double rzent(int sample, double eta) const;
  double rzext(int sample, double eta) const;
  double rmid(int sample, double eta) const;
  double rent(int sample, double eta) const;
  double rext(int sample, double eta) const;
  double zmid(int sample, double eta) const;
  double zent(int sample, double eta) const;
  double zext(int sample, double eta) const;
  double rpos(int sample,
              double eta,
              int subpos = CaloSubPos::SUBPOS_MID) const;
  double zpos(int sample,
              double eta,
              int subpos = CaloSubPos::SUBPOS_MID) const;
  double rzpos(int sample,
               double eta,
               int subpos = CaloSubPos::SUBPOS_MID) const;
  // cppcheck-suppress negativeContainerIndex; false positive
  bool isCaloBarrel(int sample) const { return m_isCaloBarrel[sample]; };
  static std::string SamplingName(int sample);

  TGraphErrors* GetGraph(unsigned int sample) const
  {
    return m_graph_layers[sample];
  };
  void SetDoGraphs(bool dographs = true) { m_dographs = dographs; };
  bool DoGraphs() const { return m_dographs; };

  TGraph* DrawGeoSampleForPhi0(int sample, int calocol, bool print = false);
  TCanvas* DrawGeoForPhi0();

  FCAL_ChannelMap* GetFCAL_ChannelMap() { return &m_FCal_ChannelMap; }
  void SetFCal_ChannelMap(const FCAL_ChannelMap* fcal_ChannnelMap)
  {
    m_FCal_ChannelMap = *fcal_ChannnelMap;
  }
  void calculateFCalRminRmax();
  virtual bool checkFCalGeometryConsistency();
  virtual void PrintMapInfo(int i, int j);

#ifdef USE_GPU
  // needed by the FCS-GPU
  t_cellmap* get_cells() { return &m_cells; };
  size_t get_n_regions(int sample)
  {
    return m_cells_in_regions[sample].size();
  };
  size_t get_region_size(int sample, int region)
  {
    return m_cells_in_regions[sample][region]->size();
  };
  int get_tot_regions()
  {
    int n = 0;
    for (int s = 0; s < MAX_SAMPLING; ++s)
      n += m_cells_in_regions[s].size();
    return n;
  }
  CaloGeometryLookup* get_region(int sample, int region)
  {
    return m_cells_in_regions[sample][region];
  };
#endif

protected:
  virtual void addcell(const CaloDetDescrElement* cell);

  virtual void post_process(int layer);

  virtual void InitRZmaps();

  t_cellmap m_cells;
  std::vector<t_cellmap> m_cells_in_sampling;
  std::vector<t_eta_cellmap> m_cells_in_sampling_for_phi0;
  std::vector<std::vector<CaloGeometryLookup*> > m_cells_in_regions;

  std::vector<bool> m_isCaloBarrel;
  std::vector<double> m_min_eta_sample[2];  //[side][calosample]
  std::vector<double> m_max_eta_sample[2];  //[side][calosample]
  std::vector<FSmap<double, double> > m_rmid_map[2];  //[side][calosample]
  std::vector<FSmap<double, double> > m_zmid_map[2];  //[side][calosample]
  std::vector<FSmap<double, double> > m_rent_map[2];  //[side][calosample]
  std::vector<FSmap<double, double> > m_zent_map[2];  //[side][calosample]
  std::vector<FSmap<double, double> > m_rext_map[2];  //[side][calosample]
  std::vector<FSmap<double, double> > m_zext_map[2];  //[side][calosample]

  bool m_dographs;
  std::vector<TGraphErrors*> m_graph_layers;
  FCAL_ChannelMap m_FCal_ChannelMap;  // for hit-to-cell assignment in FCal
  std::vector<double> m_FCal_rmin, m_FCal_rmax;

  /*
     double  m_min_eta_sample[2][MAX_SAMPLING]; //[side][calosample]
     double  m_max_eta_sample[2][MAX_SAMPLING]; //[side][calosample]
     FSmap< double , double > m_rmid_map[2][MAX_SAMPLING]; //[side][calosample]
     FSmap< double , double > m_zmid_map[2][MAX_SAMPLING]; //[side][calosample]
     FSmap< double , double > m_rent_map[2][MAX_SAMPLING]; //[side][calosample]
     FSmap< double , double > m_zent_map[2][MAX_SAMPLING]; //[side][calosample]
     FSmap< double , double > m_rext_map[2][MAX_SAMPLING]; //[side][calosample]
     FSmap< double , double > m_zext_map[2][MAX_SAMPLING]; //[side][calosample]
     */
};

#endif
