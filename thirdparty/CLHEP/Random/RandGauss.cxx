/*
  Copyright (C) 2002-2018 CERN for the benefit of the ATLAS collaboration
 */

#include "CLHEP/Random/RandGauss.h"

namespace CLHEP {

double RandGauss::shoot(HepRandomEngine *engine, double mean, double stdDev)
{
  return engine->gauss(mean, stdDev);
}

} // namespace CLHEP
