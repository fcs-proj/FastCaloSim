// Copyright (c) 2024 CERN for the benefit of the FastCaloSim project

#ifndef ISF_FASTCALOSIMEVENT_DoubleArray_h
#define ISF_FASTCALOSIMEVENT_DoubleArray_h

#include "TArrayD.h"
#include "TObject.h"

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

#endif
