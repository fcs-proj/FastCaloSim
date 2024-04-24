/*
  Copyright (C) 2002-2019 CERN for the benefit of the ATLAS collaboration
*/

#include "FastCaloSim/Geometry/CaloSampling.h"

namespace
{

const char* const sample_names[] = {
#define CALOSAMPLING(name, inbarrel, inendcap) #name,
#include "FastCaloSim/Geometry/CaloSampling.def"
#undef CALOSAMPLING
};

}  // anonymous namespace

std::string CaloSampling::getSamplingName(CaloSample theSample)
{
  return sample_names[theSample];
}

std::string CaloSampling::getSamplingName(unsigned int theSample)
{
  if (theSample >= getNumberOfSamplings())
    return "";
  return sample_names[theSample];
}

CaloSampling::CaloSample CaloSampling::getSampling(const std::string& name)
{
  for (int samp = 0; samp != static_cast<int>(Unknown); ++samp) {
    if (sample_names[samp] == name) {
      return static_cast<CaloSample>(samp);
    }
  }
  return Unknown;
}
