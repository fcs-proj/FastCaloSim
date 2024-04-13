/*
  Copyright (C) 2002-2024 CERN for the benefit of the ATLAS collaboration
*/

#include <fstream>
#include <sstream>
#include <stdexcept>

#include "FastCaloSim/VNetworkLWTNN.h"

// For reading and writing to root
#include "TFile.h"
#include "TTree.h"

VNetworkLWTNN::VNetworkLWTNN(const VNetworkLWTNN& copy_from)
    : VNetworkBase(copy_from)
    , m_json(copy_from.m_json)
    , m_printable_name(copy_from.m_printable_name)
{
  if (m_json.length() == 0) {
    throw std::invalid_argument(
        "Trying to copy a VNetworkLWTNN with length 0 m_json, probably "
        "deleteAllButNet was called on the object being copied from.");
  };
};

VNetworkLWTNN::~VNetworkLWTNN() {};

// This setup is going to do it's best to
// fill in m_json.
void VNetworkLWTNN::setupPersistedVariables()
{
  if (this->isFile(m_inputFile)) {
    ATH_MSG_DEBUG("Making an LWTNN network using a file on disk, "
                  << m_inputFile);
    m_printable_name = m_inputFile;
    fillJson();
  } else {
    ATH_MSG_DEBUG("Making an LWTNN network using a json in memory, length "
                  << m_inputFile.length());
    m_printable_name = "JSON from memory";
    m_json = m_inputFile;
  };
};

void VNetworkLWTNN::print(std::ostream& strm) const
{
  strm << m_printable_name;
};

void VNetworkLWTNN::writeNetToTTree(TTree& tree)
{
  writeStringToTTree(tree, m_json);
};

void VNetworkLWTNN::fillJson(std::string const& tree_name)
{
  ATH_MSG_VERBOSE("Trying to fill the m_json variable");
  if (this->isRootFile()) {
    ATH_MSG_VERBOSE("Treating input file as a root file");
    TFile tfile(this->m_inputFile.c_str(), "READ");
    TTree* tree = (TTree*)tfile.Get(tree_name.c_str());
    std::string found = this->readStringFromTTree(*tree);
    ATH_MSG_DEBUG("Read json from root file, length " << found.length());
    m_json = found;
  } else {
    ATH_MSG_VERBOSE("Treating input file as a text json file");
    // The input file is read into a stringstream
    std::ifstream input(m_inputFile);
    std::ostringstream sstr;
    sstr << input.rdbuf();
    m_json = sstr.str();
    input.close();
    ATH_MSG_DEBUG("Read json from text file");
  }
}

std::string VNetworkLWTNN::readStringFromTTree(TTree& tree)
{
  std::string found = std::string();
  std::string* to_found = &found;
  tree.SetBranchAddress("lwtnn_json", &to_found);
  tree.GetEntry(0);
  return found;
};

void VNetworkLWTNN::writeStringToTTree(TTree& tree, std::string json_string)
{
  tree.Branch("lwtnn_json", &json_string);
  tree.Fill();
  tree.Write();
};

void VNetworkLWTNN::deleteAllButNet()
{
  ATH_MSG_DEBUG("Replacing m_inputFile with unknown");
  m_inputFile.assign("unknown");
  m_inputFile.shrink_to_fit();
  ATH_MSG_DEBUG("Emptying the m_json string");
  m_json.clear();
  m_json.shrink_to_fit();
  ATH_MSG_VERBOSE("m_json now has capacity "
                  << m_json.capacity() << ". m_inputFile now has capacity "
                  << m_inputFile.capacity()
                  << ". m_printable_name now has capacity "
                  << m_printable_name.capacity());
};
