/*
  Copyright (C) 2002-2018 CERN for the benefit of the ATLAS collaboration
 */

#ifndef RandPoisson_h
#define RandPoisson_h 1

#include "CLHEP/Random/RandomEngine.h"

namespace CLHEP {
    
namespace RandPoisson {
  double shoot(HepRandomEngine *engine, double mean);
}

} // namespace CLHEP

#endif // RandPoisson_h
