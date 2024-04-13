/*
  Copyright (C) 2002-2018 CERN for the benefit of the ATLAS collaboration
 */

#ifndef RandFlat_h
#define RandFlat_h 1

#include "CLHEP/Random/RandomEngine.h"

namespace CLHEP {
    
namespace RandFlat {
  double shoot(HepRandomEngine *engine);
  double shoot(HepRandomEngine *engine, double a, double b);
  double shoot(HepRandomEngine *engine, double width);
}

} // namespace CLHEP

#endif // RandFlat_h
