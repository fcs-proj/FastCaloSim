/*
  Copyright (C) 2002-2020 CERN for the benefit of the ATLAS collaboration
 */

#ifndef RandGaussZiggurat_h
#define RandGaussZiggurat_h 1

#include "CLHEP/Random/RandomEngine.h"

namespace CLHEP {
    
namespace RandGaussZiggurat {
  double shoot(HepRandomEngine *engine, double mean, double stdDev);
}

} // namespace CLHEP

#endif // RandGaussZiggurat_h
