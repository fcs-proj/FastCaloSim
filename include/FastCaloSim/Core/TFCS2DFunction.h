// Copyright (c) 2024 CERN for the benefit of the FastCaloSim project

#ifndef ISF_FASTCALOSIMEVENT_TFCS2DFunction_h
#define ISF_FASTCALOSIMEVENT_TFCS2DFunction_h

#include <vector>

#include "FastCaloSim/Core/TFCSFunction.h"

class TH2;

class TFCS2DFunction : public TFCSFunction
{
public:
  TFCS2DFunction() {};
  ~TFCS2DFunction() {};

  virtual int ndim() const { return 2; };

  virtual void rnd_to_fct(float& valuex,
                          float& valuey,
                          float rnd0,
                          float rnd1) const = 0;
  virtual void rnd_to_fct(float value[], const float rnd[]) const;

  static double CheckAndIntegrate2DHistogram(const TH2* hist,
                                             std::vector<double>& integral_vec,
                                             int& first,
                                             int& last);

private:
  ClassDef(TFCS2DFunction, 1)  // TFCS2DFunction
};

#endif
