// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#ifndef ISF_FASTCALOSIMEVENT_TFCS1DFunction_h
#define ISF_FASTCALOSIMEVENT_TFCS1DFunction_h

// STL includes
#include <vector>

#include "FastCaloSim/Core/TFCSFunction.h"

class TH1;

class TFCS1DFunction : public TFCSFunction
{
public:
  TFCS1DFunction() = default;
  ~TFCS1DFunction() override = default;

  auto ndim() const -> int override { return 1; };

  /// Function gets array of random numbers rnd[] in the range [0,1) as
  /// arguments and returns function value in array value. For a n-dimensional
  /// function, value and rnd should both have n elements.
  void rnd_to_fct(float value[], const float rnd[]) const override;

  /// Function gets random number rnd in the range [0,1) as argument and returns
  /// function value
  virtual auto rnd_to_fct(double rnd) const -> double = 0;

  /// The == operator compares the content of instances.
  /// The implementation in the base class only returns true for a comparison
  /// with itself
  virtual auto operator==(const TFCS1DFunction& ref) const -> bool
  {
    return this == &ref;
  };

  static auto get_maxdev(TH1*, TH1*) -> double;

  static auto CheckAndIntegrate1DHistogram(const TH1* hist,
                                           std::vector<double>& integral_vec,
                                           int& first,
                                           int& last) -> double;

  static auto generate_histogram_random_slope(int nbinsx = 50,
                                              double xmin = 0,
                                              double xmax = 1,
                                              double zerothreshold = 0.1)
      -> TH1*;
  static auto generate_histogram_random_gauss(int nbinsx = 50,
                                              int ntoy = 100000,
                                              double xmin = 1,
                                              double xmax = 5,
                                              double xpeak = 1.5,
                                              double sigma = 0.6) -> TH1*;

private:
  ClassDef(TFCS1DFunction, 2)  // TFCS1DFunction
};

#endif
