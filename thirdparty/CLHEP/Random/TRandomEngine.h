/*
  Copyright (C) 2002-2018 CERN for the benefit of the ATLAS collaboration
 */

#ifndef TRandomEngine_h
#define TRandomEngine_h 1

#include "CLHEP/Random/RandomEngine.h"

class TRandom3;

namespace CLHEP {
    
class TRandomEngine : public HepRandomEngine {
public:
  TRandomEngine();
  virtual ~TRandomEngine();

  virtual void setSeed(long seed, int dummy=0) final;

  virtual double random() final;
  virtual double gauss(double mean, double stdDev) final;
  virtual double poisson(double mean) final;

private:
  TRandom3* m_random;
};

} // namespace CLHEP

#endif // TRandomEngine_h
