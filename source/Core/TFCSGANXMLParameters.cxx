// Copyright (c) 2024 CERN for the benefit of the FastCaloSim project

///////////////////////////////////////////////////////////////////
// TFCSGANXMLParameters.cxx, (c) ATLAS Detector software             //
///////////////////////////////////////////////////////////////////

// class header include
#include <cmath>
#include <iostream>
#include <memory>
#include <sstream>

#include "FastCaloSim/Core/TFCSGANXMLParameters.h"

#include <libxml/parser.h>

#include "TMath.h"

namespace
{
// Read an XML attribute and free the memory allocated by xmlGetProp
std::string getXmlAttr(xmlNodePtr node, const char* name)
{
  xmlChar* prop = xmlGetProp(node, BAD_CAST name);
  if (!prop)
    return {};
  std::string value(reinterpret_cast<const char*>(prop));
  xmlFree(prop);
  return value;
}

int getXmlAttrInt(xmlNodePtr node, const char* name)
{
  const std::string attribute = getXmlAttr(node, name);
  return attribute.empty() ? 0 : std::stoi(attribute);
}
}  // namespace

TFCSGANXMLParameters::TFCSGANXMLParameters() = default;

TFCSGANXMLParameters::~TFCSGANXMLParameters() = default;

void TFCSGANXMLParameters::InitialiseFromXML(
    int pid, int etaMid, const std::string& FastCaloGANInputFolderName)
{
  m_fastCaloGANInputFolderName = FastCaloGANInputFolderName;
  std::string xmlFullFileName = FastCaloGANInputFolderName + "/binning.xml";

  // RAII ownership: the document must also be freed when the numeric
  // attribute conversions below throw
  const std::unique_ptr<xmlDoc, void (*)(xmlDocPtr)> doc(
      xmlParseFile(xmlFullFileName.c_str()), xmlFreeDoc);
  if (!doc) {
    FCS_MSG_WARNING("Failed to parse XML file: " << xmlFullFileName);
    return;
  }

  for (xmlNodePtr nodeRoot = doc->children; nodeRoot != nullptr;
       nodeRoot = nodeRoot->next)
  {
    if (xmlStrEqual(nodeRoot->name, BAD_CAST "Bins")) {
      for (xmlNodePtr nodeParticle = nodeRoot->children;
           nodeParticle != nullptr;
           nodeParticle = nodeParticle->next)
      {
        if (xmlStrEqual(nodeParticle->name, BAD_CAST "Particle")) {
          int nodePid = getXmlAttrInt(nodeParticle, "pid");

          if (nodePid == pid) {
            for (xmlNodePtr nodeBin = nodeParticle->children;
                 nodeBin != nullptr;
                 nodeBin = nodeBin->next)
            {
              if (xmlStrEqual(nodeBin->name, BAD_CAST "Bin")) {
                int nodeEtaMin = getXmlAttrInt(nodeBin, "etaMin");
                int nodeEtaMax = getXmlAttrInt(nodeBin, "etaMax");
                int regionId = getXmlAttrInt(nodeBin, "regionId");

                if (std::abs(etaMid) > nodeEtaMin
                    && std::abs(etaMid) < nodeEtaMax)
                {
                  m_symmetrisedAlpha =
                      ReadBooleanAttribute("symmetriseAlpha", nodeParticle);
                  m_ganVersion = getXmlAttrInt(nodeBin, "ganVersion");
                  m_latentDim = getXmlAttrInt(nodeParticle, "latentDim");

                  for (xmlNodePtr nodeLayer = nodeBin->children;
                       nodeLayer != nullptr;
                       nodeLayer = nodeLayer->next)
                  {
                    if (xmlStrEqual(nodeLayer->name, BAD_CAST "Layer")) {
                      std::vector<double> edges;
                      std::string s(getXmlAttr(nodeLayer, "r_edges"));

                      std::istringstream ss(s);
                      std::string token;

                      while (std::getline(ss, token, ',')) {
                        if (token.empty())
                          continue;
                        edges.push_back(std::stod(token));
                      }

                      int binsInAlpha = getXmlAttrInt(nodeLayer, "n_bin_alpha");
                      int layer = getXmlAttrInt(nodeLayer, "id");

                      std::string name = "hist_pid_" + std::to_string(nodePid)
                          + "_region_" + std::to_string(regionId) + "_layer_"
                          + std::to_string(layer);
                      int xBins = static_cast<int>(edges.size()) - 1;
                      if (xBins <= 0) {
                        FCS_MSG_DEBUG(
                            "No bins defined in r for layer "
                            << layer
                            << ", setting to 1 bin to avoid empty histogram");
                        xBins = 1;
                        if (edges.empty())
                          edges.push_back(0);
                        edges.push_back(edges.back() + 1);
                      } else {
                        m_relevantlayers.push_back(layer);
                      }
                      double minAlpha = -TMath::Pi();
                      if (m_symmetrisedAlpha && binsInAlpha > 1) {
                        minAlpha = 0;
                      }
                      auto& h = m_binning[layer] = TH2D(name.c_str(),
                                                        name.c_str(),
                                                        xBins,
                                                        edges.data(),
                                                        binsInAlpha,
                                                        minAlpha,
                                                        TMath::Pi());
                      // Detach from any ROOT directory: the histogram is
                      // owned by this object and must not be registered in
                      // the (thread-unsafe) global directory
                      h.SetDirectory(nullptr);
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
}

bool TFCSGANXMLParameters::ReadBooleanAttribute(const std::string& name,
                                                xmlNodePtr node)
{
  return getXmlAttr(node, name.c_str()) == "true";
}

void TFCSGANXMLParameters::Print() const
{
  FCS_MSG_INFO("Parameters taken from XML");
  FCS_MSG_INFO("  symmetrisedAlpha: " << m_symmetrisedAlpha);
  FCS_MSG_INFO("  ganVersion:" << m_ganVersion);
  FCS_MSG_INFO("  latentDim: " << m_latentDim);
  FCS_MSG(INFO) << "  relevantlayers: ";
  for (const auto& l : m_relevantlayers) {
    FCS_MSG(INFO) << l << " ";
  }
  FCS_MSG(INFO) << END_FCS_MSG(INFO);

  for (const auto& element : m_binning) {
    int layer = element.first;
    const TH2D* h = &element.second;

    if (h->IsZombie()) {
      FCS_MSG_WARNING("Histogram for layer " << layer
                                             << " is broken. Skipping.");
      continue;
    }

    int xBinNum = h->GetNbinsX();
    const TAxis* x = h->GetXaxis();

    // If only one bin in r means layer is empty, no value should be added
    if (xBinNum == 1) {
      FCS_MSG_INFO("layer " << layer << " not used");
      continue;
    }
    FCS_MSG_INFO("Binning along r for layer " << layer);
    FCS_MSG(INFO) << "0,";
    // First fill energies
    for (int ix = 1; ix <= xBinNum; ++ix) {
      FCS_MSG(INFO) << x->GetBinUpEdge(ix) << ",";
    }
    FCS_MSG(INFO) << END_FCS_MSG(INFO);
  }
}
