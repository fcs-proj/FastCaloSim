// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

//////////////////////////////////////////////////////////////////
// TFCSGANDetectorRegion.h, (c) ATLAS Detector software
///////////////////////////////////////////////////////////////////

#pragma once

#include <map>
#include <string>
#include <vector>

// XML reader
#include <libxml/tree.h>

#include "FastCaloSim/Core/MLogging.h"
#include "TH2D.h"

class TFCSGANXMLParameters : public FastCaloSim::MLogging
{
public:
  using Binning = std::map<int, TH2D>;

  TFCSGANXMLParameters();
  ~TFCSGANXMLParameters() override;

  void InitialiseFromXML(int pid,
                         int etaMid,
                         const std::string& FastCaloGANInputFolderName);
  void Print() const;

  auto GetRelevantLayers() const -> const std::vector<int>&
  {
    return m_relevantlayers;
  };
  auto GetBinning() const -> const Binning& { return m_binning; };
  auto GetLatentSpaceSize() const -> int { return m_latentDim; };
  auto GetGANVersion() const -> int { return m_ganVersion; };
  auto IsSymmetrisedAlpha() const -> bool { return m_symmetrisedAlpha; };
  auto GetInputFolder() const -> const std::string&
  {
    return m_fastCaloGANInputFolderName;
  };

private:
  static auto ReadBooleanAttribute(const std::string& name, xmlNodePtr node)
      -> bool;

  bool m_symmetrisedAlpha {};
  Binning m_binning;
  std::vector<int> m_relevantlayers;
  int m_ganVersion {};
  int m_latentDim {};
  std::string m_fastCaloGANInputFolderName;

  ClassDefOverride(TFCSGANXMLParameters, 3)  // TFCSGANXMLParameters
};
