/*
  Copyright (C) 2002-2018 CERN for the benefit of the ATLAS collaboration
 */

#include "CLHEP/Random/RandFlat.h"

namespace CLHEP {

double RandFlat::shoot(HepRandomEngine *engine)
{
  return engine->random();
}

double RandFlat::shoot(HepRandomEngine *engine, double a, double b)
{
  return (b - a) * engine->random() + a;
}

double RandFlat::shoot(HepRandomEngine *engine, double width)
{
  return width * engine->random();
}

} // namespace CLHEP
