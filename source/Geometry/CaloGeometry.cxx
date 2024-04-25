/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#include <fstream>
#include <sstream>

#include "FastCaloSim/Geometry/CaloGeometry.h"

#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TRandom.h>
#include <TTree.h>
#include <TVector2.h>
#include <TVector3.h>

#include "FastCaloSim/Geometry/CaloDetDescrElement.h"
// #include "ISF_FastCaloSimParametrization/CaloDetDescrElement.h"
#include "FastCaloSim/Geometry/CaloSampling.h"
#include "FastCaloSim/Geometry/FastCaloSim_CaloCell_ID.h"
// #include "TMVA/Tools.h"
// #include "TMVA/Factory.h"

using namespace std;

const int CaloGeometry_calocol[24] = {
    1,   2,   3, 4,  // LAr barrel
    1,   2,   3, 4,  // LAr EM endcap
    1,   2,   3, 4,  // Hadronic end cap cal.
    1,   2,   3,  // Tile barrel
    -42, -28, 6,  // Tile gap (ITC & scint)
    1,   2,   3,  // Tile extended barrel
    1,   2,   3  // Forward EM endcap
};

const int CaloGeometry::MAX_SAMPLING =
    CaloCell_ID_FCS::MaxSample;  // number of calorimeter layers/samplings

const Identifier CaloGeometry::m_debug_identify;
std::atomic<bool> CaloGeometry::m_debug = false;

CaloGeometry::CaloGeometry()
    : m_cells_in_sampling(MAX_SAMPLING)
    , m_cells_in_sampling_for_phi0(MAX_SAMPLING)
    , m_cells_in_regions(MAX_SAMPLING)
    , m_isCaloBarrel(MAX_SAMPLING)
    , m_dographs(false)
    , m_FCal_ChannelMap(0)
{
  // TMVA::Tools::Instance();
  for (int i = 0; i < 2; ++i) {
    m_min_eta_sample[i].resize(MAX_SAMPLING);  //[side][calosample]
    m_max_eta_sample[i].resize(MAX_SAMPLING);  //[side][calosample]
    m_rmid_map[i].resize(MAX_SAMPLING);  //[side][calosample]
    m_zmid_map[i].resize(MAX_SAMPLING);  //[side][calosample]
    m_rent_map[i].resize(MAX_SAMPLING);  //[side][calosample]
    m_zent_map[i].resize(MAX_SAMPLING);  //[side][calosample]
    m_rext_map[i].resize(MAX_SAMPLING);  //[side][calosample]
    m_zext_map[i].resize(MAX_SAMPLING);  //[side][calosample]
  }
  m_graph_layers.resize(MAX_SAMPLING);
  //  for(unsigned int
  //  i=CaloCell_ID_FCS::FirstSample;i<CaloCell_ID_FCS::MaxSample;++i) {
  for (unsigned int i = CaloSampling::PreSamplerB; i <= CaloSampling::FCAL2;
       ++i)
  {
    m_graph_layers[i] = nullptr;
    CaloSampling::CaloSample s = static_cast<CaloSampling::CaloSample>(i);
    m_isCaloBarrel[i] =
        (CaloSampling::barrelPattern() & CaloSampling::getSamplingPattern(s))
        != 0;
  }
  m_isCaloBarrel[CaloCell_ID_FCS::TileGap3] = false;
}

CaloGeometry::~CaloGeometry() = default;

void CaloGeometry::addcell(const CaloDetDescrElement* cell)
{
  int sampling = cell->getSampling();
  Identifier identify = cell->identify();

  m_cells[identify] = cell;
  m_cells_in_sampling[sampling][identify] = cell;

  // m_cells[identify]= new CaloDetDescrElement(*cell);
  // m_cells_in_sampling[sampling][identify]= new CaloDetDescrElement(*cell);

  CaloGeometryLookup* lookup = nullptr;
  for (unsigned int i = 0; i < m_cells_in_regions[sampling].size(); ++i) {
    if (m_cells_in_regions[sampling][i]->IsCompatible(cell)) {
      lookup = m_cells_in_regions[sampling][i];
      // cout<<sampling<<":  found compatible map"<<endl;
      break;
    }
  }
  if (!lookup) {
    lookup = new CaloGeometryLookup(m_cells_in_regions[sampling].size());
    m_cells_in_regions[sampling].push_back(lookup);
  }
  lookup->add(cell);
}

void CaloGeometry::PrintMapInfo(int i, int j)
{
  cout << "  map " << j << ": " << m_cells_in_regions[i][j]->size() << " cells";
  if (i < 21) {
    cout << ", " << m_cells_in_regions[i][j]->cell_grid_eta() << "*"
         << m_cells_in_regions[i][j]->cell_grid_phi();
    cout << ", deta=" << m_cells_in_regions[i][j]->deta().mean() << "+-"
         << m_cells_in_regions[i][j]->deta().rms();
    cout << ", dphi=" << m_cells_in_regions[i][j]->dphi().mean() << "+-"
         << m_cells_in_regions[i][j]->dphi().rms();
    cout << ", mineta=" << m_cells_in_regions[i][j]->mineta()
         << ", maxeta=" << m_cells_in_regions[i][j]->maxeta();
    cout << ", minphi=" << m_cells_in_regions[i][j]->minphi()
         << ", maxphi=" << m_cells_in_regions[i][j]->maxphi();
    cout << endl << "         ";
    cout << ", mineta_raw=" << m_cells_in_regions[i][j]->mineta_raw()
         << ", maxeta_raw=" << m_cells_in_regions[i][j]->maxeta_raw();
    cout << ", minphi_raw=" << m_cells_in_regions[i][j]->minphi_raw()
         << ", maxphi_raw=" << m_cells_in_regions[i][j]->maxphi_raw();
    cout << endl;
  } else {
    cout << ", " << m_cells_in_regions[i][j]->cell_grid_eta() << "*"
         << m_cells_in_regions[i][j]->cell_grid_phi();
    cout << ", dx=" << m_cells_in_regions[i][j]->dx().mean() << "+-"
         << m_cells_in_regions[i][j]->dx().rms();
    cout << ", dy=" << m_cells_in_regions[i][j]->dy().mean() << "+-"
         << m_cells_in_regions[i][j]->dy().rms();
    cout << ", mindx=" << m_cells_in_regions[i][j]->mindx();
    cout << ", mindy=" << m_cells_in_regions[i][j]->mindy();
    cout << ", minx=" << m_cells_in_regions[i][j]->minx()
         << ", maxx=" << m_cells_in_regions[i][j]->maxx();
    cout << ", miny=" << m_cells_in_regions[i][j]->miny()
         << ", maxy=" << m_cells_in_regions[i][j]->maxy();
    cout << endl << "         ";
    cout << ", minx_raw=" << m_cells_in_regions[i][j]->minx_raw()
         << ", maxx_raw=" << m_cells_in_regions[i][j]->maxx_raw();
    cout << ", miny_raw=" << m_cells_in_regions[i][j]->miny_raw()
         << ", maxy_raw=" << m_cells_in_regions[i][j]->maxy_raw();
    cout << endl;
  }
}

void CaloGeometry::post_process(int sampling)
{
  // cout<<"post processing sampling "<<sampling<<endl;
  bool found_overlap = false;
  for (unsigned int j = 0; j < m_cells_in_regions[sampling].size(); ++j) {
    /*
    cout<<"Sample "<<sampling<<": checking map
    "<<j<<"/"<<m_cells_in_regions[sampling].size(); for(unsigned int
    k=0;k<m_cells_in_regions[sampling].size();++k) { cout<<",
    "<<k<<":"<<m_cells_in_regions[sampling][k]->size();
    }
    cout<<endl;
    */
    for (unsigned int i = j + 1; i < m_cells_in_regions[sampling].size(); ++i) {
      if (m_cells_in_regions[sampling][i]->has_overlap(
              m_cells_in_regions[sampling][j]))
      {
        if (!found_overlap) {
          cout << "Sample " << sampling
               << ", starting maps : " << m_cells_in_regions[sampling].size();
          for (unsigned int k = 0; k < m_cells_in_regions[sampling].size(); ++k)
          {
            cout << ", " << k << ":" << m_cells_in_regions[sampling][k]->size();
          }
          cout << endl;
        }
        found_overlap = true;
        /*
        cout<<"Sample "<<sampling<<": Found overlap between map "<<j<<" and
        "<<i<<", "
            <<m_cells_in_regions[sampling].size()<<" total maps"<<endl;
        cout<<"  current configuration map
        "<<j<<"/"<<m_cells_in_regions[sampling].size(); for(unsigned int
        k=0;k<m_cells_in_regions[sampling].size();++k) { cout<<",
        "<<k<<":"<<m_cells_in_regions[sampling][k]->size();
        }
        cout<<endl;

        PrintMapInfo(sampling,j);
        PrintMapInfo(sampling,i);
        */

        CaloGeometryLookup* toremove = m_cells_in_regions[sampling][i];
        toremove->merge_into_ref(m_cells_in_regions[sampling][j]);

        /*
        cout<<"  New ";
        PrintMapInfo(sampling,j);
        */

        for (unsigned int k = i; k < m_cells_in_regions[sampling].size() - 1;
             ++k)
        {
          m_cells_in_regions[sampling][k] = m_cells_in_regions[sampling][k + 1];
          m_cells_in_regions[sampling][k]->set_index(k);
        }
        m_cells_in_regions[sampling].resize(m_cells_in_regions[sampling].size()
                                            - 1);

        /*
        cout<<"  new configuration map
        "<<j<<"/"<<m_cells_in_regions[sampling].size(); for(unsigned int
        k=0;k<m_cells_in_regions[sampling].size();++k) { cout<<",
        "<<k<<":"<<m_cells_in_regions[sampling][k]->size();
        }
        cout<<endl;
        */

        --i;
      }
    }
  }

  if (found_overlap) {
    cout << "Sample " << sampling
         << ", final maps : " << m_cells_in_regions[sampling].size();
    for (unsigned int k = 0; k < m_cells_in_regions[sampling].size(); ++k) {
      cout << ", " << k << ":" << m_cells_in_regions[sampling][k]->size();
    }
    cout << endl;
  }

  if (found_overlap) {
    cout << "Run another round of ";
    post_process(sampling);
  }
}

void CaloGeometry::InitRZmaps()
{
  FSmap<double, double> rz_map_eta[2][MAX_SAMPLING];
  FSmap<double, double> rz_map_rmid[2][MAX_SAMPLING];
  FSmap<double, double> rz_map_zmid[2][MAX_SAMPLING];
  FSmap<double, double> rz_map_rent[2][MAX_SAMPLING];
  FSmap<double, double> rz_map_zent[2][MAX_SAMPLING];
  FSmap<double, double> rz_map_rext[2][MAX_SAMPLING];
  FSmap<double, double> rz_map_zext[2][MAX_SAMPLING];
  FSmap<double, int> rz_map_n[2][MAX_SAMPLING];

  for (unsigned int side = 0; side <= 1; ++side)
    for (unsigned int sample = 0; sample < MAX_SAMPLING; ++sample) {
      m_min_eta_sample[side][sample] = +1000;
      m_max_eta_sample[side][sample] = -1000;
    }

  for (t_cellmap::iterator calo_iter = m_cells.begin();
       calo_iter != m_cells.end();
       ++calo_iter)
  {
    const CaloDetDescrElement* theDDE = (*calo_iter).second;
    if (theDDE) {
      unsigned int sample = theDDE->getSampling();

      int side = 0;
      int sign_side = -1;
      double eta_raw = theDDE->eta_raw();
      if (eta_raw > 0) {
        side = 1;
        sign_side = +1;
      }

      if (!m_cells_in_sampling_for_phi0[sample][eta_raw]) {
        m_cells_in_sampling_for_phi0[sample][eta_raw] = theDDE;
      } else {
        if (TMath::Abs(theDDE->phi())
            < TMath::Abs(m_cells_in_sampling_for_phi0[sample][eta_raw]->phi()))
        {
          m_cells_in_sampling_for_phi0[sample][eta_raw] = theDDE;
        }
      }

      double min_eta = theDDE->eta() - theDDE->deta() / 2;
      double max_eta = theDDE->eta() + theDDE->deta() / 2;
      if (min_eta < m_min_eta_sample[side][sample])
        m_min_eta_sample[side][sample] = min_eta;
      if (max_eta > m_max_eta_sample[side][sample])
        m_max_eta_sample[side][sample] = max_eta;

      if (rz_map_eta[side][sample].find(eta_raw)
          == rz_map_eta[side][sample].end())
      {
        rz_map_eta[side][sample][eta_raw] = 0;
        rz_map_rmid[side][sample][eta_raw] = 0;
        rz_map_zmid[side][sample][eta_raw] = 0;
        rz_map_rent[side][sample][eta_raw] = 0;
        rz_map_zent[side][sample][eta_raw] = 0;
        rz_map_rext[side][sample][eta_raw] = 0;
        rz_map_zext[side][sample][eta_raw] = 0;
        rz_map_n[side][sample][eta_raw] = 0;
      }
      rz_map_eta[side][sample][eta_raw] += theDDE->eta();
      rz_map_rmid[side][sample][eta_raw] += theDDE->r();
      rz_map_zmid[side][sample][eta_raw] += theDDE->z();
      double drh = theDDE->dr() / 2;
      double dzh = theDDE->dz();
      if (sample <= CaloSampling::EMB3) {  // ensure we have a valid sampling
        drh = theDDE->dr();
      }
      // An `else` would be good here since we can't continue with a *bad*
      // sampling... Should most likely handle the situation nicely rather than
      // with a crash.

      rz_map_rent[side][sample][eta_raw] += theDDE->r() - drh;
      rz_map_zent[side][sample][eta_raw] += theDDE->z() - dzh * sign_side;
      rz_map_rext[side][sample][eta_raw] += theDDE->r() + drh;
      rz_map_zext[side][sample][eta_raw] += theDDE->z() + dzh * sign_side;
      rz_map_n[side][sample][eta_raw]++;
    }
  }

  for (int side = 0; side <= 1; ++side) {
    for (int sample = 0; sample < MAX_SAMPLING; ++sample) {
      if (!rz_map_n[side][sample].empty()) {
        for (FSmap<double, int>::iterator iter = rz_map_n[side][sample].begin();
             iter != rz_map_n[side][sample].end();
             ++iter)
        {
          double eta_raw = iter->first;
          if (iter->second < 1) {
            // ATH_MSG_WARNING("rz-map for side="<<side<<" sample="<<sample<<"
            // eta_raw="<<eta_raw<<" : #cells="<<iter->second<<" !!!");
          } else {
            double eta = rz_map_eta[side][sample][eta_raw] / iter->second;
            double rmid = rz_map_rmid[side][sample][eta_raw] / iter->second;
            double zmid = rz_map_zmid[side][sample][eta_raw] / iter->second;
            double rent = rz_map_rent[side][sample][eta_raw] / iter->second;
            double zent = rz_map_zent[side][sample][eta_raw] / iter->second;
            double rext = rz_map_rext[side][sample][eta_raw] / iter->second;
            double zext = rz_map_zext[side][sample][eta_raw] / iter->second;

            m_rmid_map[side][sample][eta] = rmid;
            m_zmid_map[side][sample][eta] = zmid;
            m_rent_map[side][sample][eta] = rent;
            m_zent_map[side][sample][eta] = zent;
            m_rext_map[side][sample][eta] = rext;
            m_zext_map[side][sample][eta] = zext;
          }
        }
        // ATH_MSG_DEBUG("rz-map for side="<<side<<" sample="<<sample<<"
        // #etas="<<m_rmid_map[side][sample].size());
      } else {
        std::cout << "rz-map for side=" << side << " sample=" << sample
                  << " is empty!!!" << std::endl;
      }
    }  // sample
  }  // side

  if (DoGraphs()) {
    for (int sample = 0; sample < MAX_SAMPLING; ++sample) {
      m_graph_layers[sample] = new TGraphErrors(rz_map_n[0][sample].size()
                                                + rz_map_n[1][sample].size());
      m_graph_layers[sample]->SetMarkerColor(
          TMath::Abs(CaloGeometry_calocol[sample]));
      m_graph_layers[sample]->SetLineColor(
          TMath::Abs(CaloGeometry_calocol[sample]));
      int np = 0;
      for (int side = 0; side <= 1; ++side) {
        for (FSmap<double, int>::iterator iter = rz_map_n[side][sample].begin();
             iter != rz_map_n[side][sample].end();
             ++iter)
        {
          double eta_raw = iter->first;
          int sign_side = -1;
          if (eta_raw > 0)
            sign_side = +1;
          // double eta =rz_map_eta[side][sample][eta_raw]/iter->second;
          double rmid = rz_map_rmid[side][sample][eta_raw] / iter->second;
          double zmid = rz_map_zmid[side][sample][eta_raw] / iter->second;
          // double rent=rz_map_rent[side][sample][eta_raw]/iter->second;
          // double zent=rz_map_zent[side][sample][eta_raw]/iter->second;
          double rext = rz_map_rext[side][sample][eta_raw] / iter->second;
          double zext = rz_map_zext[side][sample][eta_raw] / iter->second;
          m_graph_layers[sample]->SetPoint(np, zmid, rmid);
          /*
          if(isCaloBarrel(sample)) {
            m_graph_layers[sample]->SetPointError(np,0,rext-rmid);
          } else {
            m_graph_layers[sample]->SetPointError(np,(zext-zent)*sign_side,0);
          }
          */
          m_graph_layers[sample]->SetPointError(
              np, (zext - zmid) * sign_side, rext - rmid);
          ++np;
        }
      }
    }
  }
}

TGraph* CaloGeometry::DrawGeoSampleForPhi0(int sample, int calocol, bool print)
{
  if (sample < CaloSampling::PreSamplerB || sample >= CaloSampling::Unknown) {
    return nullptr;
  }
  TGraph* firstgr = nullptr;
  cout << "Start sample " << sample << " (" << SamplingName(sample) << ")"
       << endl;
  int ngr = 0;
  for (t_eta_cellmap::iterator calo_iter =
           m_cells_in_sampling_for_phi0[sample].begin();
       calo_iter != m_cells_in_sampling_for_phi0[sample].end();
       ++calo_iter)
  {
    const CaloDetDescrElement* theDDE = (*calo_iter).second;
    if (theDDE) {
      TVector3 cv;
      TGraph* gr = new TGraph(5);
      gr->SetLineColor(TMath::Abs(calocol));
      gr->SetFillColor(TMath::Abs(calocol));
      if (calocol < 0) {
        gr->SetFillStyle(1001);
      } else {
        gr->SetFillStyle(0);
      }
      gr->SetLineWidth(2);
      double r = theDDE->r();
      double dr = theDDE->dr();
      double x = theDDE->x();
      double dx = theDDE->dx();
      double y = theDDE->y();
      double dy = theDDE->dy();
      double z = theDDE->z();
      double dz = theDDE->dz() * 2;
      double eta = theDDE->eta();
      double deta = theDDE->deta();

      if (CaloSampling::PreSamplerB <= sample && sample <= CaloSampling::EMB3) {
        dr *= 2;
      }
      if (print) {
        cout << "sample=" << sample << " r=" << r << " dr=" << dr
             << " eta=" << eta << " deta=" << deta << " x=" << x << " y=" << y
             << " z=" << z << " dz=" << dz << endl;
      }
      if (isCaloBarrel(sample)) {
        cv.SetPtEtaPhi(r - dr / 2, eta - deta / 2, 0);
        gr->SetPoint(0, cv.Z(), cv.Pt());
        gr->SetPoint(4, cv.Z(), cv.Pt());
        cv.SetPtEtaPhi(r - dr / 2, eta + deta / 2, 0);
        gr->SetPoint(1, cv.Z(), cv.Pt());
        cv.SetPtEtaPhi(r + dr / 2, eta + deta / 2, 0);
        gr->SetPoint(2, cv.Z(), cv.Pt());
        cv.SetPtEtaPhi(r + dr / 2, eta - deta / 2, 0);
        gr->SetPoint(3, cv.Z(), cv.Pt());
      } else {
        if (sample < CaloSampling::FCAL0) {
          cv.SetPtEtaPhi(1, eta - deta / 2, 0);
          cv *= (z - dz / 2) / cv.Z();
          gr->SetPoint(0, cv.Z(), cv.Pt());
          gr->SetPoint(4, cv.Z(), cv.Pt());
          cv.SetPtEtaPhi(1, eta + deta / 2, 0);
          cv *= (z - dz / 2) / cv.Z();
          gr->SetPoint(1, cv.Z(), cv.Pt());
          cv.SetPtEtaPhi(1, eta + deta / 2, 0);
          cv *= (z + dz / 2) / cv.Z();
          gr->SetPoint(2, cv.Z(), cv.Pt());
          cv.SetPtEtaPhi(1, eta - deta / 2, 0);
          cv *= (z + dz / 2) / cv.Z();
          gr->SetPoint(3, cv.Z(), cv.Pt());
        } else {
          double minr = r;
          double maxr = r;
          for (double px = x - dx / 2; px <= x + dx / 2; px += dx) {
            for (double py = y - dy / 2; py <= y + dy / 2; py += dy) {
              double pr = TMath::Sqrt(px * px + py * py);
              minr = TMath::Min(minr, pr);
              maxr = TMath::Max(maxr, pr);
            }
          }
          cv.SetXYZ(minr, 0, z - dz / 2);
          gr->SetPoint(0, cv.Z(), cv.Pt());
          gr->SetPoint(4, cv.Z(), cv.Pt());
          cv.SetXYZ(maxr, 0, z - dz / 2);
          gr->SetPoint(1, cv.Z(), cv.Pt());
          cv.SetXYZ(maxr, 0, z + dz / 2);
          gr->SetPoint(2, cv.Z(), cv.Pt());
          cv.SetXYZ(minr, 0, z + dz / 2);
          gr->SetPoint(3, cv.Z(), cv.Pt());
        }
      }
      // if(calocol[sample]>0) gr->Draw("Lsame");
      //  else gr->Draw("LFsame");
      gr->Draw("LFsame");
      if (ngr == 0)
        firstgr = gr;
      ++ngr;
    }
  }
  cout << "Done sample " << sample << " (" << SamplingName(sample)
       << ")=" << ngr << endl;
  return firstgr;
}

TCanvas* CaloGeometry::DrawGeoForPhi0()
{
  TCanvas* c = new TCanvas("CaloGeoForPhi0", "Calo geometry for #phi~0");
  TH2D* hcalolayout =
      new TH2D("hcalolayoutPhi0",
               "Reconstruction geometry: calorimeter layout;z [mm];r [mm]",
               50,
               -7000,
               7000,
               50,
               0,
               4000);
  hcalolayout->Draw();
  hcalolayout->SetStats(0);
  hcalolayout->GetYaxis()->SetTitleOffset(1.4);

  TLegend* leg = new TLegend(0.30, 0.13, 0.70, 0.37);
  leg->SetFillStyle(0);
  leg->SetFillColor(10);
  leg->SetBorderSize(1);
  leg->SetNColumns(2);

  for (int sample = 0; sample < MAX_SAMPLING; ++sample) {
    TGraph* gr = DrawGeoSampleForPhi0(
        sample, CaloGeometry_calocol[sample], (sample == 17));
    if (gr) {
      std::string sname = Form("Sampling %2d : ", sample);
      sname += SamplingName(sample);
      leg->AddEntry(gr, sname.c_str(), "LF");
    }
  }
  leg->Draw();
  return c;
}

const CaloDetDescrElement* CaloGeometry::getDDE(Identifier identify)
{
  return m_cells[identify];
}
const CaloDetDescrElement* CaloGeometry::getDDE(int sampling,
                                                Identifier identify)
{
  return m_cells_in_sampling[sampling][identify];
}

const CaloDetDescrElement* CaloGeometry::getDDE(
    int sampling, float eta, float phi, float* distance, int* steps)
{
  if (sampling < 0)
    return nullptr;
  if (sampling >= MAX_SAMPLING)
    return nullptr;
  if (m_cells_in_regions[sampling].empty())
    return nullptr;

  float dist = 0;
  const CaloDetDescrElement* bestDDE = nullptr;
  if (!distance)
    distance = &dist;
  *distance = +10000000;
  int intsteps;
  int beststeps;
  if (steps)
    beststeps = (*steps);
  else
    beststeps = 0;

  if (sampling < 21) {
    for (int skip_range_check = 0; skip_range_check <= 1; ++skip_range_check) {
      for (unsigned int j = 0; j < m_cells_in_regions[sampling].size(); ++j) {
        if (!skip_range_check) {
          if (eta < m_cells_in_regions[sampling][j]->mineta())
            continue;
          if (eta > m_cells_in_regions[sampling][j]->maxeta())
            continue;
        }
        if (steps)
          intsteps = (*steps);
        else
          intsteps = 0;
        if (m_debug) {
          cout << "CaloGeometry::getDDE : check map" << j
               << " skip_range_check=" << skip_range_check << endl;
        }
        float newdist;
        const CaloDetDescrElement* newDDE =
            m_cells_in_regions[sampling][j]->getDDE(
                eta, phi, &newdist, &intsteps);
        if (m_debug) {
          cout << "CaloGeometry::getDDE : map" << j << " dist=" << newdist
               << " best dist=" << *distance << " steps=" << intsteps << endl;
        }
        if (newdist < *distance) {
          bestDDE = newDDE;
          *distance = newdist;
          if (steps)
            beststeps = intsteps;
          if (newdist < -0.1)
            break;  // stop, we are well within the hit cell
        }
      }
      if (bestDDE)
        break;
    }
  } else {
    return nullptr;
    // for(int skip_range_check=0;skip_range_check<=1;++skip_range_check) {
    // for(unsigned int j=0;j<m_cells_in_regions[sampling].size();++j) {
    // if(!skip_range_check) {
    // if(eta<m_cells_in_regions[sampling][j]->minx()) continue;
    // if(eta>m_cells_in_regions[sampling][j]->maxx()) continue;
    // if(phi<m_cells_in_regions[sampling][j]->miny()) continue;
    // if(phi>m_cells_in_regions[sampling][j]->maxy()) continue;
    //}
    // if(steps) intsteps=*steps;
    // else intsteps=0;
    // if(m_debug) {
    // cout<<"CaloGeometry::getDDE : check map"<<j<<"
    // skip_range_check="<<skip_range_check<<endl;
    //}
    // float newdist;
    // const CaloGeoDetDescrElement*
    // newDDE=m_cells_in_regions[sampling][j]->getDDE(eta,phi,&newdist,&intsteps);
    // if(m_debug) {
    // cout<<"CaloGeometry::getDDE : map"<<j<<" dist="<<newdist<<" best
    // dist="<<*distance<<" steps="<<intsteps<<endl;
    //}
    // if(newdist<*distance) {
    // bestDDE=newDDE;
    //*distance=newdist;
    // if(steps) beststeps=intsteps;
    // if(newdist<-0.1) break; //stop, we are well within the hit cell
    //}
    //}
    // if(bestDDE) break;
    //}
  }
  if (steps)
    *steps = beststeps;
  return bestDDE;
}

const CaloDetDescrElement* CaloGeometry::getFCalDDE(
    int sampling, float x, float y, float z, float* distance, int* steps)
{
  int isam = sampling - 20;
  int iphi(-100000), ieta(-100000);
  Long64_t mask1[] {0x34, 0x34, 0x35};
  Long64_t mask2[] {0x36, 0x36, 0x37};
  bool found = m_FCal_ChannelMap.getTileID(isam, x, y, ieta, iphi);
  if (steps && found)
    *steps = 0;
  if (!found) {
    // cout << "Warning: Hit is not matched with any FCal cell! Looking for the
    // closest cell" << endl;
    found = getClosestFCalCellIndex(sampling, x, y, ieta, iphi, steps);
  }
  if (!found) {
    cout << "Error: Unable to find the closest FCal cell!" << endl;
    return nullptr;
  }

  // cout << "CaloGeometry::getFCalDDE: x:" << x << " y: " << y << " ieta: " <<
  // ieta << " iphi: " << iphi << " cmap->x(): " <<
  // m_FCal_ChannelMap.x(isam,ieta,iphi) << " cmap->y(): " <<
  // m_FCal_ChannelMap.y(isam,ieta,iphi) << endl;
  Long64_t id = (ieta << 5) + 2 * iphi;
  if (isam == 2)
    id += (8 << 8);

  if (z > 0)
    id += (mask2[isam - 1] << 12);
  else
    id += (mask1[isam - 1] << 12);

  id = id << 44;
  Identifier identify((unsigned long long)id);
  const CaloDetDescrElement* foundcell = m_cells[identify];

  // Report the shortest distance in x or y to the face of the cell if outside
  // the cell. If inside the cell this will be negative and be the shortest
  // distance to a cell face.
  float dist = 0;
  if (!distance)
    distance = &dist;
  float distanceX = abs(foundcell->x() - x) - foundcell->dx() / 2;
  float distanceY = abs(foundcell->y() - y) - foundcell->dy() / 2;
  *distance = max(distanceX, distanceY);

  return foundcell;
}

bool CaloGeometry::getClosestFCalCellIndex(
    int sampling, float x, float y, int& ieta, int& iphi, int* steps)
{
  double rmin = m_FCal_rmin[sampling - 21];
  double rmax = m_FCal_rmax[sampling - 21];
  int isam = sampling - 20;
  double a = 1.;
  const double b = 0.01;
  const int nmax = 100;
  int i = 0;

  const double r = sqrt(x * x + y * y);
  if (r == 0.)
    return false;
  const double r_inverse = 1. / r;

  if ((r / rmax) > (rmin * r_inverse)) {
    x = x * rmax * r_inverse;
    y = y * rmax * r_inverse;
    while ((!m_FCal_ChannelMap.getTileID(isam, a * x, a * y, ieta, iphi))
           && i < nmax)
    {
      a -= b;
      i++;
    }
  } else {
    x = x * rmin * r_inverse;
    y = y * rmin * r_inverse;
    while ((!m_FCal_ChannelMap.getTileID(isam, a * x, a * y, ieta, iphi))
           && i < nmax)
    {
      a += b;
      i++;
    }
  }
  if (steps)
    *steps = i + 1;
  return i < nmax;
}

bool CaloGeometry::PostProcessGeometry()
{
  for (int i = 0; i < MAX_SAMPLING; ++i) {
    post_process(i);
    for (unsigned int j = 0; j < m_cells_in_regions[i].size(); ++j) {
      m_cells_in_regions[i][j]->post_process();
    }
    // if(i>=21) break;
  }

  InitRZmaps();

  /*
  cout<<"all : "<<m_cells.size()<<endl;
  for(int sampling=0;sampling<MAX_SAMPLING;++sampling) {
    cout<<"cells sampling "<<sampling<<" :
  "<<m_cells_in_sampling[sampling].size()<<" cells"; cout<<",
  "<<m_cells_in_regions[sampling].size()<<" lookup map(s)"<<endl; for(unsigned
  int j=0;j<m_cells_in_regions[sampling].size();++j) { PrintMapInfo(sampling,j);
      //break;
    }
    //if(i>0) break;
  }
  */

  return true;
}

void CaloGeometry::Validate(int nrnd)
{
  int ntest = 0;
  cout << "start CaloGeometry::Validate()" << endl;
  for (t_cellmap::iterator ic = m_cells.begin(); ic != m_cells.end(); ++ic) {
    const CaloDetDescrElement* cell = ic->second;
    int sampling = cell->getSampling();
    // if(sampling>=21) continue;

    if (m_debug_identify == cell->identify()) {
      cout << "CaloGeometry::Validate(), cell " << ntest
           << " id=" << cell->identify() << endl;
      m_debug = true;
    }

    for (int irnd = 0; irnd < nrnd; ++irnd) {
      std::stringstream pos;
      std::stringstream cellpos;
      std::stringstream foundcellpos;
      int steps = 0;
      float distance = 0;
      bool is_inside;
      bool is_inside_foundcell = false;
      const CaloDetDescrElement* foundcell = nullptr;
      if (sampling < 21) {
        float eta = cell->eta() + 1.95 * (gRandom->Rndm() - 0.5) * cell->deta();
        float phi = cell->phi() + 1.95 * (gRandom->Rndm() - 0.5) * cell->dphi();
        foundcell = getDDE(sampling, eta, phi, &distance, &steps);

        pos << "eta=" << eta << " phi=" << phi;
        cellpos << "eta=" << cell->eta() << " eta_raw=" << cell->eta_raw()
                << " deta=" << cell->deta() << " ("
                << (cell->eta_raw() - cell->eta()) / cell->deta() << ") ; "
                << "phi=" << cell->phi() << " phi_raw=" << cell->phi_raw()
                << " dphi=" << cell->dphi() << " ("
                << (cell->phi_raw() - cell->phi()) / cell->dphi() << ")";
        if (foundcell) {
          foundcellpos << "eta=" << foundcell->eta()
                       << " eta_raw=" << foundcell->eta_raw()
                       << " deta=" << foundcell->deta() << " ("
                       << (foundcell->eta_raw() - foundcell->eta())
                  / foundcell->deta()
                       << ") ; "
                       << "phi=" << foundcell->phi()
                       << " phi_raw=" << foundcell->phi_raw()
                       << " dphi=" << foundcell->dphi() << " ("
                       << (foundcell->phi_raw() - foundcell->phi())
                  / cell->dphi()
                       << ")";
          is_inside_foundcell =
              TMath::Abs((eta - foundcell->eta()) / foundcell->deta()) < 0.55
              && TMath::Abs((phi - foundcell->phi()) / foundcell->dphi())
                  < 0.55;
        }
        is_inside = TMath::Abs((eta - cell->eta()) / cell->deta()) < 0.49
            && TMath::Abs((phi - cell->phi()) / cell->dphi()) < 0.49;
      } else {
        float x = cell->x() + 1.95 * (gRandom->Rndm() - 0.5) * cell->dx();
        float y = cell->y() + 1.95 * (gRandom->Rndm() - 0.5) * cell->dy();
        float z = cell->z();
        foundcell = getFCalDDE(sampling, x, y, z);

        pos << "x=" << x << " y=" << y << " z=" << z;
        cellpos << "x=" << cell->x() << " x_raw=" << cell->x_raw()
                << " dx=" << cell->dx() << " ("
                << (cell->x_raw() - cell->x()) / cell->dx() << ") ; "
                << "y=" << cell->y() << " y_raw=" << cell->y_raw()
                << " dy=" << cell->dy() << " ("
                << (cell->y_raw() - cell->y()) / cell->dy() << ")";
        if (foundcell) {
          foundcellpos << "x=" << foundcell->x()
                       << " x_raw=" << foundcell->x_raw()
                       << " dx=" << foundcell->dx() << " ("
                       << (foundcell->x_raw() - foundcell->x())
                  / foundcell->dx()
                       << ") ; "
                       << "y=" << foundcell->y()
                       << " y_raw=" << foundcell->y_raw()
                       << " dy=" << foundcell->dy() << " ("
                       << (foundcell->y_raw() - foundcell->y()) / cell->dy()
                       << ")";
          is_inside_foundcell =
              TMath::Abs((x - foundcell->x()) / foundcell->dx()) < 0.75
              && TMath::Abs((y - foundcell->y()) / foundcell->dy()) < 0.75;
        }
        is_inside = TMath::Abs((x - cell->x()) / cell->dx()) < 0.49
            && TMath::Abs((y - cell->y()) / cell->dy()) < 0.49;
        // m_debug=true;
      }

      if (m_debug && foundcell) {
        cout << "CaloGeometry::Validate(), irnd=" << irnd
             << ": cell id=" << cell->identify() << ", sampling=" << sampling
             << ", foundcell id=" << foundcell->identify() << ", " << steps
             << " steps" << endl;
        cout << "  " << cellpos.str() << endl;
      }
      if (cell == foundcell) {
        if (steps > 3 && distance < -0.01) {
          cout << "cell id=" << cell->identify() << ", sampling=" << sampling
               << " found in " << steps << " steps, dist=" << distance << " "
               << pos.str() << endl;
          cout << "  " << cellpos.str() << endl;
        }
      } else {
        if (!foundcell) {
          cout << "cell id=" << cell->identify() << " not found!";
          cout << " No cell found in " << steps << " steps, dist=" << distance
               << " " << pos.str() << endl;
          cout << "  input sampling=" << sampling << " " << cellpos.str()
               << endl;
          // return;
        }
        if (is_inside && foundcell && !is_inside_foundcell) {
          cout << "cell id=" << cell->identify()
               << " not found, but inside cell area!";
          cout << " Found instead id=" << foundcell->identify() << " in "
               << steps << " steps, dist=" << distance << " " << pos.str()
               << endl;
          cout << "  input sampling=" << sampling << " " << cellpos.str()
               << endl;
          cout << "  output sampling=" << foundcell->getSampling() << " "
               << foundcellpos.str() << endl;
          return;
        }
      }
    }
    m_debug = false;
    if (ntest % 25000 == 0)
      cout << "Validate cell " << ntest << " with " << nrnd << " random hits"
           << endl;
    ++ntest;
  }
  cout << "end CaloGeometry::Validate()" << endl;
}

double CaloGeometry::deta(int sample, double eta) const
{
  int side = 0;
  if (eta > 0)
    side = 1;

  double mineta = m_min_eta_sample[side][sample];
  double maxeta = m_max_eta_sample[side][sample];

  if (eta < mineta) {
    return fabs(eta - mineta);
  } else if (eta > maxeta) {
    return fabs(eta - maxeta);
  } else {
    double d1 = fabs(eta - mineta);
    double d2 = fabs(eta - maxeta);
    if (d1 < d2)
      return -d1;
    else
      return -d2;
  }
}

void CaloGeometry::minmaxeta(int sample,
                             double eta,
                             double& mineta,
                             double& maxeta) const
{
  int side = 0;
  if (eta > 0)
    side = 1;

  mineta = m_min_eta_sample[side][sample];
  maxeta = m_max_eta_sample[side][sample];
}

double CaloGeometry::rmid(int sample, double eta) const
{
  int side = 0;
  if (eta > 0)
    side = 1;

  return m_rmid_map[side][sample].find_closest(eta)->second;
}

double CaloGeometry::zmid(int sample, double eta) const
{
  int side = 0;
  if (eta > 0)
    side = 1;

  return m_zmid_map[side][sample].find_closest(eta)->second;
}

double CaloGeometry::rzmid(int sample, double eta) const
{
  int side = 0;
  if (eta > 0)
    side = 1;

  if (isCaloBarrel(sample))
    return m_rmid_map[side][sample].find_closest(eta)->second;
  else
    return m_zmid_map[side][sample].find_closest(eta)->second;
}

double CaloGeometry::rent(int sample, double eta) const
{
  int side = 0;
  if (eta > 0)
    side = 1;

  return m_rent_map[side][sample].find_closest(eta)->second;
}

double CaloGeometry::zent(int sample, double eta) const
{
  int side = 0;
  if (eta > 0)
    side = 1;

  return m_zent_map[side][sample].find_closest(eta)->second;
}

double CaloGeometry::rzent(int sample, double eta) const
{
  int side = 0;
  if (eta > 0)
    side = 1;

  if (isCaloBarrel(sample))
    return m_rent_map[side][sample].find_closest(eta)->second;
  else
    return m_zent_map[side][sample].find_closest(eta)->second;
}

double CaloGeometry::rext(int sample, double eta) const
{
  int side = 0;
  if (eta > 0)
    side = 1;

  return m_rext_map[side][sample].find_closest(eta)->second;
}

double CaloGeometry::zext(int sample, double eta) const
{
  int side = 0;
  if (eta > 0)
    side = 1;

  return m_zext_map[side][sample].find_closest(eta)->second;
}

double CaloGeometry::rzext(int sample, double eta) const
{
  int side = 0;
  if (eta > 0)
    side = 1;

  if (isCaloBarrel(sample))
    return m_rext_map[side][sample].find_closest(eta)->second;
  else
    return m_zext_map[side][sample].find_closest(eta)->second;
}

double CaloGeometry::rpos(int sample, double eta, int subpos) const
{
  int side = 0;
  if (eta > 0)
    side = 1;

  if (subpos == CaloSubPos::SUBPOS_ENT)
    return m_rent_map[side][sample].find_closest(eta)->second;
  if (subpos == CaloSubPos::SUBPOS_EXT)
    return m_rext_map[side][sample].find_closest(eta)->second;

  return m_rmid_map[side][sample].find_closest(eta)->second;
}

double CaloGeometry::zpos(int sample, double eta, int subpos) const
{
  int side = 0;
  if (eta > 0)
    side = 1;

  if (subpos == CaloSubPos::SUBPOS_ENT)
    return m_zent_map[side][sample].find_closest(eta)->second;
  if (subpos == CaloSubPos::SUBPOS_EXT)
    return m_zext_map[side][sample].find_closest(eta)->second;
  return m_zmid_map[side][sample].find_closest(eta)->second;
}

double CaloGeometry::rzpos(int sample, double eta, int subpos) const
{
  int side = 0;
  if (eta > 0)
    side = 1;

  if (isCaloBarrel(sample)) {
    if (subpos == CaloSubPos::SUBPOS_ENT)
      return m_rent_map[side][sample].find_closest(eta)->second;
    if (subpos == CaloSubPos::SUBPOS_EXT)
      return m_rext_map[side][sample].find_closest(eta)->second;
    return m_rmid_map[side][sample].find_closest(eta)->second;
  } else {
    if (subpos == CaloSubPos::SUBPOS_ENT)
      return m_zent_map[side][sample].find_closest(eta)->second;
    if (subpos == CaloSubPos::SUBPOS_EXT)
      return m_zext_map[side][sample].find_closest(eta)->second;
    return m_zmid_map[side][sample].find_closest(eta)->second;
  }
}

std::string CaloGeometry::SamplingName(int sample)
{
  return CaloSampling::getSamplingName(sample);
}

void CaloGeometry::calculateFCalRminRmax()
{
  m_FCal_rmin.resize(3, FLT_MAX);
  m_FCal_rmax.resize(3, 0.);

  double x(0.), y(0.), r(0.);
  for (int imap = 1; imap <= 3; imap++)
    for (auto it = m_FCal_ChannelMap.begin(imap);
         it != m_FCal_ChannelMap.end(imap);
         it++)
    {
      x = it->second.x();
      y = it->second.y();
      r = sqrt(x * x + y * y);
      if (r < m_FCal_rmin[imap - 1])
        m_FCal_rmin[imap - 1] = r;
      if (r > m_FCal_rmax[imap - 1])
        m_FCal_rmax[imap - 1] = r;
    }
}

bool CaloGeometry::checkFCalGeometryConsistency()
{
  unsigned long long phi_index, eta_index;
  float x, y, dx, dy;
  long id;

  long mask1[] {0x34, 0x34, 0x35};
  long mask2[] {0x36, 0x36, 0x37};

  m_FCal_rmin.resize(3, FLT_MAX);
  m_FCal_rmax.resize(3, 0.);

  for (int imap = 1; imap <= 3; imap++) {
    int sampling = imap + 20;

    if ((int)m_cells_in_sampling[sampling].size()
        != 2
            * std::distance(m_FCal_ChannelMap.begin(imap),
                            m_FCal_ChannelMap.end(imap)))
    {
      cout << "Error: Incompatibility between FCalChannel map and GEO file: "
              "Different number of cells in m_cells_in_sampling and "
              "FCal_ChannelMap"
           << endl;
      cout << "m_cells_in_sampling: " << m_cells_in_sampling[sampling].size()
           << endl;
      cout << "FCal_ChannelMap: "
           << 2
              * std::distance(m_FCal_ChannelMap.begin(imap),
                              m_FCal_ChannelMap.end(imap))
           << endl;
      return false;
    }

    for (auto it = m_FCal_ChannelMap.begin(imap);
         it != m_FCal_ChannelMap.end(imap);
         it++)
    {
      phi_index = it->first & 0xffff;
      eta_index = it->first >> 16;
      x = it->second.x();
      y = it->second.y();
      m_FCal_ChannelMap.tileSize(imap, eta_index, phi_index, dx, dy);

      id = (mask1[imap - 1] << 12) + (eta_index << 5) + 2 * phi_index;

      if (imap == 2)
        id += (8 << 8);

      Identifier id1((unsigned long long)(id << 44));
      const CaloDetDescrElement* DDE1 = getDDE(id1);

      id = (mask2[imap - 1] << 12) + (eta_index << 5) + 2 * phi_index;
      if (imap == 2)
        id += (8 << 8);
      Identifier id2((unsigned long long)(id << 44));
      const CaloDetDescrElement* DDE2 = getDDE(id2);

      if (!TMath::AreEqualRel(x, DDE1->x(), 1.E-8)
          || !TMath::AreEqualRel(y, DDE1->y(), 1.E-8)
          || !TMath::AreEqualRel(x, DDE2->x(), 1.E-8)
          || !TMath::AreEqualRel(y, DDE2->y(), 1.E-8))
      {
        cout << "Error: Incompatibility between FCalChannel map and GEO file \n"
             << x << " " << DDE1->x() << " " << DDE2->x() << y << " "
             << DDE1->y() << " " << DDE2->y() << endl;
        return false;
      }
    }
  }

  return true;
}
