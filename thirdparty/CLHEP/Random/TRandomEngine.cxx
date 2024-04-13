/*
  Copyright (C) 2002-2018 CERN for the benefit of the ATLAS collaboration
 */

#include <TRandom3.h>

#include "CLHEP/Random/TRandomEngine.h"

namespace CLHEP {

TRandomEngine::TRandomEngine()
  : HepRandomEngine()
{
  m_random = new TRandom3(42);
}

TRandomEngine::~TRandomEngine()
{
  delete m_random;
}

void TRandomEngine::setSeed(long seed, int)
{
  m_random->SetSeed(seed);
}

double TRandomEngine::random()
{
  // The use of 1 - engine->Rndm() is a fudge for TRandom3, as it generates random numbers
  // in (0,1], but [0,1) or (0,1) is needed.
  return 1 - m_random->Rndm();
}

double TRandomEngine::gauss(double mean, double stdDev)
{
  return m_random->Gaus(mean, stdDev);
}

double TRandomEngine::poisson(double mean)
{
  return m_random->Poisson(mean);
}

} // namespace CLHEP
