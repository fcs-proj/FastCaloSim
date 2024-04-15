/*
  Copyright (C) 2002-2018 CERN for the benefit of the ATLAS collaboration
 */

#ifndef RandGauss_h
#define RandGauss_h 1

#include "CLHEP/Random/RandomEngine.h"

namespace CLHEP {
    
namespace RandGauss {
  double shoot(HepRandomEngine *engine, double mean, double stdDev);
}

} // namespace CLHEP

#endif // RandGauss_h
