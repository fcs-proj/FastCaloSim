// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#ifndef ISF_FASTCALOSIMEVENT_IntArray_h
#define ISF_FASTCALOSIMEVENT_IntArray_h

#include "TArrayI.h"
#include "TObject.h"

class IntArray
    : public TObject
    , public TArrayI
{
public:
  IntArray();
  IntArray(int);
  ~IntArray() override;

private:
  ClassDef(IntArray, 1)
};

#endif
