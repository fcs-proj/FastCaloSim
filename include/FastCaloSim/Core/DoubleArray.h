// Copyright (c) 2024 CERN for the benefit of the FastCaloSim project

#ifndef DoubleArray_h
#define DoubleArray_h

#include "TArrayD.h"
#include "TObject.h"

namespace FastCaloSim::Core
{

class DoubleArray
    : public TObject
    , public TArrayD
{
public:
  DoubleArray();
  DoubleArray(int);
  ~DoubleArray();

private:
  ClassDef(DoubleArray, 1)
};

}  // namespace FastCaloSim::Core

#endif
