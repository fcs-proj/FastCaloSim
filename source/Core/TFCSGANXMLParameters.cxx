// Copyright (c) 2024 CERN for the benefit of the FastCaloSim project

///////////////////////////////////////////////////////////////////
// TFCSGANXMLParameters.cxx, (c) ATLAS Detector software             //
///////////////////////////////////////////////////////////////////

// class header include
#include <iostream>
#include <sstream>

#include "FastCaloSim/Core/TFCSGANXMLParameters.h"

#include "TMath.h"

using namespace FastCaloSim::Core;

TFCSGANXMLParameters::TFCSGANXMLParameters() {}

TFCSGANXMLParameters::~TFCSGANXMLParameters() {}

void TFCSGANXMLParameters::InitialiseFromXML(
    int pid, int etaMid, const std::string& FastCaloGANInputFolderName)
{
  m_fastCaloGANInputFolderName = FastCaloGANInputFolderName;
  std::string xmlFullFileName = FastCaloGANInputFolderName + "/binning.xml";

  xmlDocPtr doc = xmlParseFile(xmlFullFileName.c_str());
  for (xmlNodePtr nodeRoot = doc->children; nodeRoot != nullptr;
       nodeRoot = nodeRoot->next)
  {
    if (xmlStrEqual(nodeRoot->name, BAD_CAST "Bins")) {
      for (xmlNodePtr nodeParticle = nodeRoot->children;
           nodeParticle != nullptr;
           nodeParticle = nodeParticle->next)
      {
        if (xmlStrEqual(nodeParticle->name, BAD_CAST "Particle")) {
          int nodePid =
              atof((const char*)xmlGetProp(nodeParticle, BAD_CAST "pid"));
          for (xmlNodePtr nodeBin = nodeParticle->children; nodeBin != nullptr;
               nodeBin = nodeBin->next)
          {
            if (nodePid == pid) {
              if (xmlStrEqual(nodeBin->name, BAD_CAST "Bin")) {
                int nodeEtaMin =
                    atof((const char*)xmlGetProp(nodeBin, BAD_CAST "etaMin"));
                int nodeEtaMax =
                    atof((const char*)xmlGetProp(nodeBin, BAD_CAST "etaMax"));
                int regionId =
                    atof((const char*)xmlGetProp(nodeBin, BAD_CAST "regionId"));

                if (fabs(etaMid) > nodeEtaMin && fabs(etaMid) < nodeEtaMax) {
                  m_symmetrisedAlpha =
                      ReadBooleanAttribute("symmetriseAlpha", nodeParticle);
                  m_ganVersion = atof(
                      (const char*)xmlGetProp(nodeBin, BAD_CAST "ganVersion"));
                  m_latentDim = atof((const char*)xmlGetProp(
                      nodeParticle, BAD_CAST "latentDim"));

                  for (xmlNodePtr nodeLayer = nodeBin->children;
                       nodeLayer != nullptr;
                       nodeLayer = nodeLayer->next)
                  {
                    if (xmlStrEqual(nodeLayer->name, BAD_CAST "Layer")) {
                      std::vector<double> edges;
                      std::string s((const char*)xmlGetProp(
                          nodeLayer, BAD_CAST "r_edges"));

                      std::istringstream ss(s);
                      std::string token;

                      while (std::getline(ss, token, ',')) {
                        edges.push_back(atof(token.c_str()));
                      }

                      int binsInAlpha = atof((const char*)xmlGetProp(
                          nodeLayer, BAD_CAST "n_bin_alpha"));
                      int layer = atof(
                          (const char*)xmlGetProp(nodeLayer, BAD_CAST "id"));

                      std::string name = "hist_pid_" + std::to_string(nodePid)
                          + "_region_" + std::to_string(regionId) + "_layer_"
                          + std::to_string(layer);
                      int xBins = edges.size() - 1;
                      if (xBins == 0)
                        xBins = 1;  // remove warning
                      else
                        m_relevantlayers.push_back(layer);
                      double minAlpha = -TMath::Pi();
                      if (m_symmetrisedAlpha && binsInAlpha > 1) {
                        minAlpha = 0;
                      }
                      m_binning[layer] = TH2D(name.c_str(),
                                              name.c_str(),
                                              xBins,
                                              &edges[0],
                                              binsInAlpha,
                                              minAlpha,
                                              TMath::Pi());
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  xmlFreeDoc(doc);
}

bool TFCSGANXMLParameters::ReadBooleanAttribute(const std::string& name,
                                                xmlNodePtr node)
{
  std::string attribute = (const char*)xmlGetProp(node, BAD_CAST name.c_str());
  bool value = attribute == "true" ? true : false;
  return value;
}

void TFCSGANXMLParameters::Print() const
{
  MSG_INFO("Parameters taken from XML");
  MSG_INFO("  symmetrisedAlpha: " << m_symmetrisedAlpha);
  MSG_INFO("  ganVersion:" << m_ganVersion);
  MSG_INFO("  latentDim: " << m_latentDim);
  MSG(INFO) << "  relevantlayers: ";
  for (auto l : m_relevantlayers) {
    MSG(INFO) << l << " ";
  }
  MSG(INFO) << END_MSG(INFO);

  for (auto element : m_binning) {
    int layer = element.first;
    TH2D* h = &element.second;

    int xBinNum = h->GetNbinsX();
    TAxis* x = (TAxis*)h->GetXaxis();

    // If only one bin in r means layer is empty, no value should be added
    if (xBinNum == 1) {
      MSG_INFO("layer " << layer << " not used");
      continue;
    }
    MSG_INFO("Binning along r for layer " << layer);
    MSG(INFO) << "0,";
    // First fill energies
    for (int ix = 1; ix <= xBinNum; ++ix) {
      MSG(INFO) << x->GetBinUpEdge(ix) << ",";
    }
    MSG(INFO) << END_MSG(INFO);
  }
}
