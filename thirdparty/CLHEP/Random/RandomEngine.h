/*
  Copyright (C) 2002-2018 CERN for the benefit of the ATLAS collaboration
 */

#ifndef HepRandomEngine_h
#define HepRandomEngine_h 1

namespace CLHEP {
    
class HepRandomEngine {
public:
  HepRandomEngine();
  virtual ~HepRandomEngine();

  virtual void setSeed(long seed, int extra = 0) = 0;

  virtual double random() = 0;
  virtual double gauss(double mean, double stdDev) = 0;
  virtual double poisson(double mean) = 0;
};

} // namespace CLHEP

#endif // HepRandomEngine_h
