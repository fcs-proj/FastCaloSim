/*
  Copyright (C) 2002-2018 CERN for the benefit of the ATLAS collaboration
 */

#include "CLHEP/Random/RandPoisson.h"

namespace CLHEP {

double RandPoisson::shoot(HepRandomEngine *engine, double mean)
{
  return engine->poisson(mean);
}

} // namespace CLHEP
