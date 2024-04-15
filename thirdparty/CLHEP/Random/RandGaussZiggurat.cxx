/*
  Copyright (C) 2002-2020 CERN for the benefit of the ATLAS collaboration
 */

#include "CLHEP/Random/RandGaussZiggurat.h"

namespace CLHEP {

double RandGaussZiggurat::shoot(HepRandomEngine *engine, double mean, double stdDev)
{
  return engine->gauss(mean, stdDev);
}

} // namespace CLHEP
