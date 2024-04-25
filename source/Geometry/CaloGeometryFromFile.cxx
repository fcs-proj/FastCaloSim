/*
  Copyright (C) 2002-2023 CERN for the benefit of the ATLAS collaboration
*/

#include <fstream>
#include <sstream>

#include <TAxis.h>
#include <TFile.h>
#include <TGraph.h>
#include <TTree.h>

#ifdef ENABLE_XROOTD_SUPPORT
#  include "XrdStreamBuf.h"
#endif

#include "FastCaloSim/Geometry/CaloDetDescrElement.h"
#include "FastCaloSim/Geometry/CaloGeometryFromFile.h"

CaloGeometryFromFile::CaloGeometryFromFile()
    : CaloGeometry()
{
}

bool CaloGeometryFromFile::LoadGeometryFromFile(const std::string& fileName,
                                                const std::string& treeName,
                                                const std::string& hashFileName)
{
  std::map<uint64_t, uint64_t> cellId_vs_cellHashId_map;

  std::unique_ptr<std::istream> hashStream {};
  std::unique_ptr<std::streambuf> hashStreamBuf {};
#ifdef ENABLE_XROOTD_SUPPORT
  if (hashFileName.find("root://") != std::string::npos) {
    hashStreamBuf = std::make_unique<XrdStreamBuf>(hashFileName);
    hashStream = std::make_unique<std::istream>(hashStreamBuf.get());
  } else {
#endif
    std::unique_ptr<std::ifstream> hashStreamDirect =
        std::make_unique<std::ifstream>(hashFileName);
    if (!hashStreamDirect->is_open()) {
      std::cout << "Error: Could not open " << hashFileName << std::endl;
      throw std::runtime_error("Could not open file");
    }
    hashStream = std::move(hashStreamDirect);
#ifdef ENABLE_XROOTD_SUPPORT
  }
#endif

  std::cout << "Loading cellId_vs_cellHashId_map" << std::endl;

  int i = 0;
  uint64_t id, hash_id;
  while (!hashStream->eof()) {
    i++;

    *hashStream >> id >> hash_id;
    cellId_vs_cellHashId_map[id] = hash_id;
    if (i % 10000 == 0)
      std::cout << "Line: " << i << " id " << std::hex << id << " hash_id "
                << std::dec << hash_id << std::endl;
  }

  std::cout << "Done." << std::endl;

  auto f = std::unique_ptr<TFile>(TFile::Open(fileName.c_str()));
  if (!f) {
    std::cerr << "Error: Could not open file '" << fileName << "'" << std::endl;
    return false;
  }
  TTree* tree = (TTree*)f->Get(treeName.c_str());
  if (!tree)
    return false;

  TTree* fChain = tree;

  CaloDetDescrElement cell;

  // List of branches
  TBranch* b_identifier;  //!
  TBranch* b_calosample;  //!
  TBranch* b_eta;  //!
  TBranch* b_phi;  //!
  TBranch* b_r;  //!
  TBranch* b_eta_raw;  //!
  TBranch* b_phi_raw;  //!
  TBranch* b_r_raw;  //!
  TBranch* b_x;  //!
  TBranch* b_y;  //!
  TBranch* b_z;  //!
  TBranch* b_x_raw;  //!
  TBranch* b_y_raw;  //!
  TBranch* b_z_raw;  //!
  TBranch* b_deta;  //!
  TBranch* b_dphi;  //!
  TBranch* b_dr;  //!
  TBranch* b_dx;  //!
  TBranch* b_dy;  //!
  TBranch* b_dz;  //!

  fChain->SetMakeClass(1);
  fChain->SetBranchAddress("identifier", &cell.m_identify, &b_identifier);
  fChain->SetBranchAddress("calosample", &cell.m_calosample, &b_calosample);
  fChain->SetBranchAddress("eta", &cell.m_eta, &b_eta);
  fChain->SetBranchAddress("phi", &cell.m_phi, &b_phi);
  fChain->SetBranchAddress("r", &cell.m_r, &b_r);
  fChain->SetBranchAddress("eta_raw", &cell.m_eta_raw, &b_eta_raw);
  fChain->SetBranchAddress("phi_raw", &cell.m_phi_raw, &b_phi_raw);
  fChain->SetBranchAddress("r_raw", &cell.m_r_raw, &b_r_raw);
  fChain->SetBranchAddress("x", &cell.m_x, &b_x);
  fChain->SetBranchAddress("y", &cell.m_y, &b_y);
  fChain->SetBranchAddress("z", &cell.m_z, &b_z);
  fChain->SetBranchAddress("x_raw", &cell.m_x_raw, &b_x_raw);
  fChain->SetBranchAddress("y_raw", &cell.m_y_raw, &b_y_raw);
  fChain->SetBranchAddress("z_raw", &cell.m_z_raw, &b_z_raw);
  fChain->SetBranchAddress("deta", &cell.m_deta, &b_deta);
  fChain->SetBranchAddress("dphi", &cell.m_dphi, &b_dphi);
  fChain->SetBranchAddress("dr", &cell.m_dr, &b_dr);
  fChain->SetBranchAddress("dx", &cell.m_dx, &b_dx);
  fChain->SetBranchAddress("dy", &cell.m_dy, &b_dy);
  fChain->SetBranchAddress("dz", &cell.m_dz, &b_dz);

  Long64_t nentries = fChain->GetEntriesFast();
  for (Long64_t jentry = 0; jentry < nentries; jentry++) {
    Long64_t ientry = fChain->LoadTree(jentry);
    if (ientry < 0)
      break;
    fChain->GetEntry(jentry);

    if (cellId_vs_cellHashId_map.find(cell.m_identify)
        != cellId_vs_cellHashId_map.end())
    {
      cell.m_hash_id = cellId_vs_cellHashId_map[cell.m_identify];
      if (cell.m_hash_id != jentry)
        std::cout << jentry << " : ERROR hash=" << cell.m_hash_id << std::endl;
    } else
      std::cout << std::endl
                << "ERROR: Cell id not found in the cellId_vs_cellHashId_map!!!"
                << std::endl
                << std::endl;

    const CaloDetDescrElement* pcell = new CaloDetDescrElement(cell);
    this->addcell(pcell);

    if (jentry % 25000 == 0) {
      std::cout << "Checking loading cells from file" << std::endl
                << jentry << " : " << pcell->getSampling() << ", "
                << pcell->identify() << std::endl;
    }
  }

  f->Close();
  bool ok = PostProcessGeometry();

  std::cout << "Result of PostProcessGeometry(): " << ok << std::endl;

  const CaloDetDescrElement* mcell = 0;
  unsigned long long cellid64(3179554531063103488);
  Identifier cellid(cellid64);
  mcell = this->getDDE(cellid);  // This is working also for the FCal

  std::cout << "\n \n";
  std::cout << "Testing whether CaloGeoGeometry is loaded properly"
            << std::endl;
  if (!mcell) {
    std::cout << "Cell is not found" << std::endl;
  } else {
    std::cout << "Identifier " << mcell->identify() << " sampling "
              << mcell->getSampling() << " eta: " << mcell->eta()
              << " phi: " << mcell->phi() << " CaloDetDescrElement=" << mcell
              << std::endl
              << std::endl;

    const CaloDetDescrElement* mcell2 =
        this->getDDE(mcell->getSampling(), mcell->eta(), mcell->phi());
    std::cout << "Identifier " << mcell2->identify() << " sampling "
              << mcell2->getSampling() << " eta: " << mcell2->eta()
              << " phi: " << mcell2->phi() << " CaloDetDescrElement=" << mcell2
              << std::endl
              << std::endl;
  }

  return ok;
}

bool CaloGeometryFromFile::LoadFCalGeometryFromFiles(
    const std::array<std::string, 3>& fileNames)
{
  std::vector<std::unique_ptr<std::istream>> electrodes;
  std::vector<std::unique_ptr<std::streambuf>> electrodesBuf;
  electrodes.reserve(3);
  electrodesBuf.reserve(3);

  for (uint16_t i = 0; i < 3; i++) {
    const std::string& file = fileNames[i];
#ifdef ENABLE_XROOTD_SUPPORT
    if (file.find("root://") != std::string::npos) {
      electrodesBuf.emplace_back(std::make_unique<XrdStreamBuf>(file));
      electrodes.emplace_back(
          std::make_unique<std::istream>(electrodesBuf.back().get()));
    } else {
#endif
      std::unique_ptr<std::ifstream> directStream =
          std::make_unique<std::ifstream>(file);
      if (!directStream->is_open()) {
        std::cout << "Error: Could not open " << file << std::endl;
        throw std::runtime_error("Could not open file");
      }
      electrodes.push_back(std::move(directStream));
#ifdef ENABLE_XROOTD_SUPPORT
    }
#endif
  }

  int thisTubeId;
  int thisTubeI;
  int thisTubeJ;
  // int    thisTubeID;
  // int    thisTubeMod;
  double thisTubeX;
  double thisTubeY;
  TString tubeName;

  // int second_column;
  std::string seventh_column;
  std::string eight_column;
  int ninth_column;

  int i;
  for (int imodule = 1; imodule <= 3; imodule++) {
    std::cout << "Loading FCal electrode #" << imodule << std::endl;

    i = 0;
    while (1) {
      (*electrodes[imodule - 1]) >> tubeName;
      if (electrodes[imodule - 1]->eof())
        break;
      (*electrodes[imodule - 1]) >> thisTubeId;  // ?????
      (*electrodes[imodule - 1]) >> thisTubeI;
      (*electrodes[imodule - 1]) >> thisTubeJ;
      (*electrodes[imodule - 1]) >> thisTubeX;
      (*electrodes[imodule - 1]) >> thisTubeY;
      (*electrodes[imodule - 1]) >> seventh_column;
      (*electrodes[imodule - 1]) >> eight_column;
      (*electrodes[imodule - 1]) >> ninth_column;

      tubeName.ReplaceAll("'", "");
      std::string tubeNamestring = tubeName.Data();

      std::istringstream tileStream1(std::string(tubeNamestring, 1, 1));
      std::istringstream tileStream2(std::string(tubeNamestring, 3, 2));
      std::istringstream tileStream3(std::string(tubeNamestring, 6, 3));
      int a1 = 0, a2 = 0, a3 = 0;
      if (tileStream1)
        tileStream1 >> a1;
      if (tileStream2)
        tileStream2 >> a2;
      if (tileStream3)
        tileStream3 >> a3;

      std::stringstream s;

      m_FCal_ChannelMap.add_tube(tubeNamestring,
                                 imodule,
                                 thisTubeId,
                                 thisTubeI,
                                 thisTubeJ,
                                 thisTubeX,
                                 thisTubeY,
                                 seventh_column);

      i++;
    }
  }

  m_FCal_ChannelMap.finish();  // Creates maps

  electrodes.clear();
  electrodesBuf.clear();

  this->calculateFCalRminRmax();
  return this->checkFCalGeometryConsistency();
}

void CaloGeometryFromFile::DrawFCalGraph(int isam, int color)
{
  std::stringstream ss;
  ss << "FCal" << isam - 20 << std::endl;

  const int size = m_cells_in_sampling[isam].size();

  std::vector<double> x;
  std::vector<double> y;
  x.reserve(size);
  y.reserve(size);

  for (auto it = m_cells_in_sampling[isam].begin();
       it != m_cells_in_sampling[isam].end();
       it++)
  {
    x.push_back(it->second->x());
    y.push_back(it->second->y());
  }

  TGraph* graph = new TGraph(size, &x[0], &y[0]);
  graph->SetLineColor(color);
  graph->SetTitle(ss.str().c_str());
  graph->SetMarkerStyle(21);
  graph->SetMarkerColor(color);
  graph->SetMarkerSize(0.5);
  graph->GetXaxis()->SetTitle("x [mm]");
  graph->GetYaxis()->SetTitle("y [mm]");

  graph->Draw("AP");
}

void CaloGeometryFromFile::calculateFCalRminRmax()
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

bool CaloGeometryFromFile::checkFCalGeometryConsistency()
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
      std::cout
          << "Error: Incompatibility between FCalChannel map and GEO file: "
             "Different number of cells in m_cells_in_sampling and "
             "FCal_ChannelMap"
          << std::endl;
      std::cout << "m_cells_in_sampling: "
                << m_cells_in_sampling[sampling].size() << std::endl;
      std::cout << "FCal_ChannelMap: "
                << 2
              * std::distance(m_FCal_ChannelMap.begin(imap),
                              m_FCal_ChannelMap.end(imap))
                << std::endl;
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
        std::cout
            << "Error: Incompatibility between FCalChannel map and GEO file \n"
            << x << " " << DDE1->x() << " " << DDE2->x() << y << " "
            << DDE1->y() << " " << DDE2->y() << std::endl;
        return false;
      }
    }
  }

  return true;
}
