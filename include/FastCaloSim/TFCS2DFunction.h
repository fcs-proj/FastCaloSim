/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#ifndef ISF_FASTCALOSIMEVENT_TFCS2DFunction_h
#define ISF_FASTCALOSIMEVENT_TFCS2DFunction_h

#include "FastCaloSim/TFCSFunction.h"
#include <vector>

class TH2;

class TFCS2DFunction : public TFCSFunction {
public:
  TFCS2DFunction(){};
  ~TFCS2DFunction(){};

  virtual int ndim() const { return 2; };

  virtual void rnd_to_fct(float &valuex, float &valuey, float rnd0,
                          float rnd1) const = 0;
  virtual void rnd_to_fct(float value[], const float rnd[]) const;

  static double CheckAndIntegrate2DHistogram(const TH2 *hist,
                                             std::vector<double> &integral_vec,
                                             int &first, int &last);

  static void unit_test(
      TH2 *hist = nullptr, TFCS2DFunction *rtof = nullptr,
      const char *outfilename = "TFCS2DFunction_unit_test.root",
      int nrnd = 10000000);
  static void unit_tests(
      TH2 *hist = nullptr,
      const char *outfilename = "TFCS2DFunction_unit_test.root",
      int nrnd = 10000000);

private:
  ClassDef(TFCS2DFunction, 1) // TFCS2DFunction
};

#endif
