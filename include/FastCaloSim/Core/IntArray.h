// Copyright (c) 2024 CERN for the benefit of the FastCaloSim project

#ifndef IntArray_h
#define IntArray_h

#include "TArrayI.h"
#include "TObject.h"

namespace FastCaloSim::Core
{
class IntArray
    : public TObject
    , public TArrayI
{
public:
  IntArray();
  IntArray(int);
  ~IntArray();

private:
  ClassDef(IntArray, 1)
};
}  // namespace FastCaloSim::Core

#endif
