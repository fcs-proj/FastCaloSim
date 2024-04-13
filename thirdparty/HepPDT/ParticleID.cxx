#include "HepPDT/ParticleID.hh"

#include <cmath>
#include <iostream>
#include <limits>

namespace HepPDT
{

ParticleID::ParticleID(const int pdgID)
{
  if (pdgID == 11 || pdgID == 211 || pdgID == 2212)
    m_charge = 1.;
  else if (pdgID == -11 || pdgID == -211 || pdgID == -2212)
    m_charge = -1.;
  else if (pdgID == 22 || std::abs(pdgID) == 2112 || pdgID == 111)
    m_charge = 0;
  else {
    std::cerr << "Error: This pdgID is not supported: " << std::endl;
    m_charge = std::numeric_limits<double>::quiet_NaN();
  }
}

} // namespace HepPDT
